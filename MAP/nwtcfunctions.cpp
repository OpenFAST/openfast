/**
 * ====================================================================================================
 *                              NWTCFunctions.cpp
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


/**
 * ====================================================================================================
 * MSQS_Init
 *  
 * This for-loop stores each line in the MAP input file into CableLibrary,
 * Element, and Node structures. Line line character arrays are broken into
 * individual words and then converted to integers or real numbers as needed.
 * 
 * The cable library should not have parameters start with '#'; if this is
 * the case, an error is returned.
 * ====================================================================================================
 */
void 
MSQS_Init( MAP_InitInputType_class       &Init           , 
           MAP_InputType_class           &In             , 
           MAP_ParameterType_class       &Param          , 
           void*                         NULL_ContState  ,
           void*                         NULL_DiscrState ,
           MAP_ConstraintStateType_class &Constrnt       , 
           MAP_OtherStateType_class      &Data           , 
           MAP_OutputType_class          &Out            , 
           void*                         NULL_TimeInt    ,
           MAP_InitOutputType_class      &InitOut        ,
           MAP_ErrStat                   &Error          , 
           MAP_Message                   &Msg )
{     
    
#ifdef NO_PYTHON
    checkpoint();
#endif /* NO_PYTHON */

  Error.clean_error( ); // set MAP error status (presuming here there is no error
                        // just by simply extering this function from the glue code)    

  Msg.messageClean( );  // clean the Msg parameter

  Data.setMessageReferenceToUserData( Msg );
  Data.setErrorStatusReferenceToUserData( Error );

  // Each string split is stored as a vector of words
  std::vector<std::string> words;

  /**
   * ==========   Populate run-time variables   ===============     <-------------------------------+
   *  -- Element                                                                        //          |
   *  -- Node                                                                           //          |
   *  -- Cable Library                                                                  //          |
   */                                                                                   //          |
  try {                                                                                 //          |
    // Set depth environmental property                                                 //          |
    // If double depth is uninitialized:                                                //          |
    if( Init.getDepth_dbl( ) == -999.9 ) {                                              //          |
      Data.addDepth( Init.getDepth( ) , Error, Msg );  // water depth                   //          |
    }                                                                                   //          |
    else {                                                                              //          |
      Data.addDepth( Init.getDepth_dbl( ) );                                            //          |
    };                                                                                  //          |
                                                                                        //          |
    // Set sea desnity environmental property                                           //          |
    // If double sea_density is uninitialized:                                          //          |
    if( Init.getSeaDensity_dbl( ) == -999.9) {                                          //          |
      Data.addSeaDensity( Init.getSeaDensity( ) , Error, Msg );  // water density       //          |
    }                                                                                   //          |
    else {                                                                              //          |
      Data.addSeaDensity( Init.getSeaDensity_dbl( ) );                                  //          |
    };                                                                                  //          |
                                                                                        //          |
    // Set gravity environmental property                                               //          |
    // If double gravity is uninitialized:                                              //          |
    if( Init.getGravity_dbl( ) == -999.9) {                                             //          |
      Data.addGravity( Init.getGravity( ), Error, Msg );  // gravity                    //          |
    }                                                                                   //          |
    else {                                                                              //          |
      Data.addGravity( Init.getGravity_dbl( ) );                                        //          |
    };                                                                                  //          |
                                                                                      //          |
                                                                                      //          |
    /**                                                                               //          |
     * ==========   Populate CableLibrary   ===============     <----------------+    //          |
     * Here we set the following properties in the                   //          |    //          |
     *                                                               //          |    //          |
     * CableLibrary structure                                        //          |    //          |
     */                                                              //          |    //          |
    for( int i=0 ; i<Init.sizeOfCableLibrary() ; i++){               //          |    //          |
      // This takes strings from the input file and separates      //          |    //          |
      // them into individual words                                //          |    //          |
      boost::split( words                                    ,     //          |    //          |  
                    Init.getCableLibraryData( i, Error, Msg ),     //          |    //          |
                    boost::is_any_of( " \n" )                ,     //          |    //          |
                    boost::token_compress_on );                    //          |    //          |
                                                                   //          |    //          |
      // make sure the number of words in the string is 6          //          |    //          |
      if ( words.size( ) != 6 ) throw MAP_ERROR_5;                 //          |    //          |
                                                                   //          |    //          |
      // No take the words and pass them to the user data class.   //          |    //          |
      //  This class will sort them into individual variables      //          |    //          |
      Data.addCableLibrary( words , i+1 , Error, Msg );            //          |    //          |
                                                                   //          |    //          |
      words.clear();                                               //          |    //          |  
    };// END for                                                     //   -------+    //          |
    //============== <END> =======================================================    //          |
    //          |
    //          |
    /**                                                                               //          |
     * ==========   Populate Node   ===============     <------------------------+    //          |
     * Here we set the following properties in the                   //          |    //          |
     * Node object:                                                  //          |    //          |
     *     X  - position in the global coordinate                    //          |    //          |
     *     Y  - position in the global coordinate                    //          |    //          |
     *     Z  - position in the global coordinate                    //          |    //          |
     *     M  - node mass                                            //          |    //          |
     *     B  - node displacement (volume)                           //          |    //          |
     *     FX - applied force in X                                   //          |    //          |
     *     FY - applied force in Y                                   //          |    //          |
     *     FZ - applied force in Z                                   //          |    //          |
     */                                                              //          |    //          |
    for( int i=0 ; i<Init.sizeOfNodeData() ; i++) {                  //          |    //          |
      // This takes strings from the input file and separates      //          |    //          |
      // them into individual words                                //          |    //          |         
      boost::split( words                   ,                      //          |    //          | 
                    Init.getNodeData( i )   ,                      //          |    //          |
                    boost::is_any_of(" \n") ,                      //          |    //          | 
                    boost::token_compress_on );                    //          |    //          | 
      //          |    //          |
      // make sure the number of words in the string is 11         //          |    //          |
      if ( words.size( ) !=11 ) throw MAP_ERROR_13;                //          |    //          |
      //          |    //          |
      // create node based on the string in Init.getNodeData( )    //          |    //          |
      Data.addNode( words , i, Error, Msg );                       //          |    //          |
      //          |    //          |
      // clear words for the next pass in the for-loop             //          |    //          |
      words.clear();                                               //          |    //          |
    };// END for                                                     //   -------+    //          |
    //============== <END> =======================================================    //          |
    //          |
    //          |
    /**                                                                               //          |
     * ==========   Populate Element   ===============     <---------------------+    //          |
     * This initilizes the Element data structure, which by and      //          |    //          |
     * large, is an object of pointers to other variables help       //          |    //          |
     * throughout the MAP program.                                   //          |    //          |
     */                                                              //          |    //          |
    //          |    //          |
    for( int i=0 ; i<Init.sizeOfElementData() ; i++) {               //          |    //          |
      //          |    //          |
      // This takes strings from the input file and separates them //          |    //          |
      // into individual words                                     //          |    //          |
      boost::split( words                    ,                     //          |    //          |
                    Init.getElementData( i ) ,                     //          |    //          |
                    boost::is_any_of(" \n")  ,                     //          |    //          |
                    boost::token_compress_on );                    //          |    //          |
      //          |    //          |
      // make sure the number of words in the string is at least 4 //          |    //          |
      if ( words.size( ) <= 4 ) throw MAP_ERROR_20;                //          |    //          |
      //          |    //          |
      Data.addElement( words, i, Error, Msg );                     //          |    //          |
      //          |    //          |
      words.clear();                                               //          |    //          |
    };// END for                                                     //   -------+    //          |
    //============== <END> =======================================================    //          |
    //          |
  } catch ( MAP_ERROR_CODE code ) {                                                     //          |
    std::string str = "";                                                             //          |
    str += boost::lexical_cast < std::string > ( code );                              //          |
    if ( str.size( ) == 1 ) str = " " + str;                                          //          |
    str += "] : " + MAP_ERROR_CODE_to_string.at( code );                              //          |
    //          |
    Error.set_error_key( MAP_ERROR );                                                 //          |
    Msg.RecordErrorToErrorList( str );                                                            //          |
    //          |
    str.erase(2,1);                                                                   //          |
    std::ostringstream S;                                                             //          |
    S.str("");S.clear();                                                              //          |
    S << ">>>> " +str +"    :    MAP error in file ";                                 //          |
    S << __FILE__;                                                                    //          |
    S << " on line ";                                                                 //          |
    S << __LINE__;                                                                    //          |
    S << " occurred.";                                                                //          |
    Msg.WriteErrorToOutputFile( S.str() );                                                       //          |
  }// END try                                                                           //   -------+
  //============== <END> ============================================================================



  /**
   * ==========   Create associations between Node and Elements   ===============     <-------------+
   *  -- Element                                                                        //          |
   *  -- Node                                                                           //          |
   *  -- Cable Library                                                                  //          |
   */                                                                                   //          |
                                                                                        //          |
  // As long as no errors are raised in the node, element and cable library             //          |
  // initialization phase, then we can proceed with the next steps in the process.      //          |
  if ( Error.error_status() == MAP_SAFE || Error.error_status() == MAP_WARNING ) {      //          |
                                                                                        //          |
    /**                                                                               //          |
     * ==========   set NWTC types  ===============     <--------------------+        //          |
     * MAP_ParameterType_class                                         //          |        //          |
     * MAP_ConstraintStateType_class                                   //          |        //          |
     * MAP_InputType_class                                             //          |        //          |
     * MAP_OutputType_class                                            //          |        //          |
     *                                                           //          |        //          |
     *        RULES                                              //          |        //          |
     * Parameter types are those that are :                      //          |        //          |
     *    - all variable in CableLibrary                         //          |        //          |
     *    - all variables in Element and NodeS                   //          |        //          |
     *      which are time-invariant (i.e.,                      //          |        //          |
     *      is_fixed = true )                                    //          |        //          |
     *                                                           //          |        //          |
     *  Constraint types are those that are :                    //          |        //          |
     *   - all variables in Element and Nodes                    //          |        //          |
     *     which are time varying (i.e.,                         //          |        //          |
     *     is_fixed = false )                                    //          |        //          |
     */                                                          //          |        //          |
    //          |        //          |
    // Assign data in the CableLibrary, assign as a              //          |        //          |
    // MAP_ParameterType_class (this CANNOT be changed)                //          |        //          |
    for( int i=0 ; i<Data.getSizeOfCableLibrary() ; i++){        //          |        //          |
      Data.setCableLibraryReference( Param, i );               //          |        //          |
    };//END for                                                  //----------+        //          |
    //============== <END> ===================================================        //          |
    //          |
    //          |
    /**                                                                               //          |
     * ==========   Node VarType assignment  ===============     <------------+       //          |
     * Assign data in the Node                                    //          |       //          |
     */                                                           //          |       //          |
    for( int i=0 ; i<Data.getSizeOfNode( ) ; i++) {               //          |       //          |
      // if the node is a 'Connect', Newton's equation          //          |       //          |
      // must be written for it                                 //          |       //          |
      if ( Data.getNodeType( i ) == Connect ) {                 //          |       //          |
        // This sets a reference to this node in              //          |       //          |
        // UserData. UserData is a structure used             //          |       //          |
        // by the Numerics class to calculate                 //          |       //          |
        // the equations we are solving                       //          |       //          |
        Data.setNodeReferenceToUserData( i );                 //          |       //          |
      };// END if                                               //          |       //          |
      //          |       //          |
      // Assign X, Y, Z, FX, FY, FZ, M and B as a parameter,    //          |       //          |
      // constraint or input type in MAP                        //          |       //          |
      Data.associateNodeVarTypeToNWTCType( i        ,           //          |       //          |
                                           In       ,           //          |       //          |
                                           Param    ,           //          |       //          |
                                           Constrnt ,           //          |       //          |
                                           Out      ,           //          |       //          |
                                           Error    ,           //          |       //          |
                                           Msg);                //          |       //          |
    };//END for                                                   //----------+       //          |
    //============== <END> ====================================================       //          |
    //          |
    //          |
    /**                                                                               //          |
     * ==========   Element VarType assignment  ============     <------------+       //          |
     * Assign data in the Element                                 //          |       //          |
     *                                                            //          |       //          |
     * (only Lu applied in this case, since line_property,        //          |       //          |
     * anchor and fairlead are pointers to variables in           //          |       //          |
     * MAP_OtherStateType_class manipulated by other types Assign       //          |       //          |
     * data in the CableLibrary, assign as a MAP_ParameterType_class    //          |       //          |
     * (default)                                                  //          |       //          |
     */                                                           //          |       //          |
    //          |       //          |
    for( int i=0 ; i<Data.getSizeOfElement() ; i++) {             //          |       //          |
      //          |       //          |
      // Each element has 2 equations associated with it        //          |       //          |
      Data.incrementNumEquations();                             //          |       //          |
      Data.incrementNumEquations();                             //          |       //          |
      //          |       //          |
      // This sets a reference to this node in                  //          |       //          |
      // UserData. UserData is a structure used                 //          |       //          |
      // by the Numerics class to calculate                     //          |       //          |
      // the equations we are solving                           //          |       //          |
      Data.setElementReferenceToUserData( i );                  //          |       //          |
      //          |       //          |
      // Assign Lu, H and V as a parameter,                     //          |       //          |
      // constraint or input type in MAP                        //          |       //          |
      Data.associateElementVarTypeToNWTCType( i        ,        //          |       //          |
                                              Param    ,        //          |       //          |
                                              Constrnt ,        //          |       //          |
                                              Error    ,        //          |       //          |
                                              Msg );            //          |       //          |
    };//END for                                                   //----------+       //          |
    //============== <END> ====================================================       //          |
                                                                                      //          |
    // set UserData reference                                                         //          |
                                                                                      //          |
    // This sets a reference to this node in UserData. UserData is a structure used   //          |
    // by the Numerics class to calculate the equations we are solving                //          |
    Data.setMAP_ConstraintStateType_classReferenceToUserData( Constrnt );                   //          |
                                                                                      //          |
    /**                                                                               //          |
     * ==========   Node & Element reference check  ============     <--------+       //          |
     * As a last step, we check the reference_counter for         //          |       //          |
     * all VarType (in Node Element and CableLibrary) to          //          |       //          |
     * ensure it is pointing to one MAP_ParameterType_class,            //          |       //          |
     * MAP_ConstraintStateType_class or MAP_InputType_class.                  //          |       //          |
     */                                                           //          |       //          |
                                                                  //          |       //          |
    for( int i=0 ; i<Data.getSizeOfNode() ; i++) {                //          |       //          |
      Data.checkNodeVarTypeReferences( i , Error, Msg );          //          |       //          |
    };//END for                                                   //          |       //          |
                                                                  //          |       //          |
    for( int i=0 ; i<Data.getSizeOfElement() ; i++) {             //          |       //          |
      Data.checkElementVarTypeReferences( i , Error, Msg );       //          |       //          |
    };//END for                                                   //----------+       //          |
    //============== <END> ====================================================       //          |
             //          |
    //          |
    /**                                                                               //          |
     * ==========   node and element initial conditions  ============     <---+       //          |
     * Once the initialization is successfully completed, set     //          |       //          |
     * the initial conditions for the elements and nodes          //          |       //          |
     */                                                           //          |       //          |
                           //          |       //          |
    // set Node.sum_FX, Node.sum_FY and Node.sum_FZ               //          |       //          |
    // to zero                                                    //          |       //          |
    Data.initializeSumForce();                                    //          |       //          |
    //          |       //          |
    // This sets the local variables for each element             //          |       //          |
    // (with respect to the element frame)H and V set FX,         //          |       //          |
    for( int i=0 ; i<Data.getSizeOfElement() ; i++) {             //          |       //          |
      // we initialize the element to set:                      //          |       //          |
      //     - the angle psi (orientation of xy relative        //          |       //          |
      //       to XYZ)                                          //          |       //          |
      //     - cable area and weight per unit length            //          |       //          |
      //     - the vertical and horizontal cable excursion      //          |       //          |
      //     - the horizontal and vertical forces applied at    //          |       //          |
      //       both the fairlead and anchor                     //          |       //          |
      Data.initializeCableElement( i , Error , Msg);            //          |       //          |
    };//END for                                                   //          |       //          |
    //          |       //          |
    // this sets the initial conditions for Fix and Vessel nodes  //          |       //          |
    // This adds sum_FX, sum_FY and sum_FZ into the node          //          |       //          |
    for( int i=0 ; i<Data.getSizeOfNode() ; i++) {                //          |       //          |
      Data.initializeCableNode( i );                            //          |       //          |
    };//END for                                                   //----------+       //          |
    //============== <END> ====================================================       //          |
    //          |
    //          |
    /**                                                                               //          |
     * ==========   Initialize numerical solver   ============     <----------+       //          |
     */                                                           //          |       //          |
    for( int i=0 ; i<Init.numOfSolverOptions( ) ; i++) {          //          |       //          |
      words.clear( );                                           //          |       //          |
      boost::split( words                      ,                //          |       //          |
                    Init.getSolverOptions( i ) ,                //          |       //          |
                    boost::is_any_of(" \n")    ,                //          |       //          |
                    boost::token_compress_on );                 //          |       //          |
      //          |       //          |
      for( unsigned int j=0 ; j<words.size()-1 ; j++ ) {        //          |       //          |
        if (words[j] != " " ) {                               //          |       //          |
          // simply passes the string on the solver         //          |       //          |
          // option portion of the MAP input file into      //          |       //          |
          // a string in MAP_OtherStateType_class. This is        //          |       //          |
          // eventually used by Numerics class to           //          |       //          |
          // separate individual options                    //          |       //          |
          Data.setSolverOptions( words[j], Msg );           //          |       //          |
        };//END if                                            //          |       //          |
      };// END if                                               //          |       //          |
      //          |       //          |
      words.clear();                                            //          |       //          |
    };//END for                                                   //----------+       //          |
    //============== <END> ====================================================       //          |
    //          |
  }// END try                                                                           //   -------+
  //============== <END> ============================================================================


  /**
   * ==========   Initialize the Numerics routine   ===============     <---------------------------+
   * This routine solves the initial conditions for the quasi-static cable.             //          |
   * This is done once. After solving, the results are printed to an output             //          |
   * file, profile.map                                                                  //          |
   */                                                                                   //          |
  //          |
  // Only initialize the solver if no errors occurred                                   //          |
  //          |
  // print statics results to a file "profile.map"                                      //          |
  // Get current date/time, format is YYYY-MM-DD.HH:mm:ss                               //          |
  time_t     now = time(0);                                                             //          |
  struct tm  tstruct;                                                                   //          |
  char       buf[80];                                                                   //          |
  tstruct = *localtime(&now);                                                           //          |
  //          |
  strftime( buf , sizeof(buf) , "%Y-%m-%d.%X", &tstruct);                               //          |
  //          |
  std::ofstream out_file("summary.map");                                                //          |
  out_file << "MAP Profile File.\n";                                                    //          |
  out_file << "Outputs were generated using MAP v1.0 on <";                             //          |
  out_file << buf;                                                                      //          |
  out_file << ">\n\n";                                                                  //          |
  //          |
  try {                                                                                 //          |
    if ( Error.error_status() == MAP_SAFE || Error.error_status() == MAP_WARNING ) {  //          |
      Data.initializeNumericSolver( Init, Error, Msg );                             //          |
                                                                                    //          |
      // call the PETSc soliver only if it has been initialized in the MAP_Init     //          |
      // routine. The PETSc numerics routine will not be initialized if the         //          |
      // '-help' flag is raised as an option in the 'NUMERIC OPTIONS' section       //          |
      // of the MAP input file                                                      //          |
                                                                                    //          |
      if ( Data.isNumericsUninitialized()==false ){                                 //          |
        Data.Solve( Error, Msg );                                                   //          |
                                                                                    //          |
        //Data.CheckResidualConvergence( Error , Msg );
        
        // this sets the initial conditions for the Fix and Vessel nodes          //          |
        // This adds sum_FX, sum_FY and sum_FZ into the node                      //          |
        for( int i=0 ; i<Data.getSizeOfNode() ; i++) {                            //          |
          Data.initializeCableNode              ( i );                          //          |
          Data.initializeEquilibriumNodePosition( i );                          //          |
        };//END for                                                               //          |
        //          |
        out_file << Data.summary();                                                  //          |
        out_file.close();                                                         //          |
      }                                                                             //          |
      else {                                                                        //          |
        throw MAP_ERROR_51;                                                       //          |
      };//END if                                                                    //          |
    }                                                                                 //          |
    else {                                                                            //          |
      throw MAP_ERROR_51;                                                           //          |
    };//END if                                                                        //          |
  } catch ( MAP_ERROR_CODE &code ) {                                                    //          |
    std::string str = "";                                                             //          |
    str += boost::lexical_cast < std::string > ( code );                              //          |
    str += "] : " + MAP_ERROR_CODE_to_string.at( code );                              //          |
    Msg.RecordErrorToErrorList( str );                                                            //          |
    //          |
    str = "";                                                                         //          |
    str += boost::lexical_cast < std::string > ( MAP_ERROR_52 );                      //          |
    str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_52 );                      //          |
    Msg.RecordErrorToErrorList( str );                                                            //          |
    //          |
    Error.set_error_key( MAP_ERROR );                                                 //          |
    //          |
    str.erase(2,1);                                                                   //          |
    std::ostringstream S;                                                             //          |
    S.str("");S.clear();                                                              //          |
    S << ">>>> " +str;
    Msg.WriteErrorToOutputFile( S.str() );                                                       //          |
    //          |
    out_file << "*************************************************************\n";    //          |
    out_file << "*                    WARNING                                *\n";    //          |
    out_file << "* The model is not initialized because the numerics routine *\n";    //          |
    out_file << "* could not be initialized (possibly due to conflicting     *\n";    //          |
    out_file << "* inputs in the MAP input file).                            *\n";    //          |
    out_file << "*************************************************************\n\n";  //          |
    out_file << Msg.GetStatusString( );                                               //          |
    out_file << "\n\n";                                                               //          |
    out_file << Data.summary();                                                         //          |
  }// END try                                                                           //   -------+
  //============== <END> ============================================================================
 
  out_file.close();

  //Data.plot(Error,Msg);
};


/**
 * ====================================================================================================
 * MSQS_UpdateStates
 * ====================================================================================================
 */
void MSQS_UpdateStates( float T                           ,
                        int interval                      ,
                        MAP_InputType_class           &In       ,
                        MAP_ParameterType_class       &Param    ,
                        void* NULL_ContinuousState        ,
                        void* NULL_DiscreteState          ,
                        MAP_ConstraintStateType_class &Constrnt ,
                        MAP_OtherStateType_class      &Data     ,
                        MAP_ErrStat             &Error    ,
                        MAP_Message             &Msg ){
    
  // clean the Msg parameter
  Msg.messageClean(); 

  // call the PETSc soliver only if it has been initialized in the MAP_Init routine.
  // The PETSc numerics routine will not be initialized if the '-help' flag is raised
  // as an option in the 'NUMERIC OPTIONS' section of the MAP input file
  if ( Data.isNumericsUninitialized()==false && Error.error_status()==MAP_SAFE ){
    // Set MAP error status (presuming here there is no error 
    // just by simply extering this function from the glue code)    
    Error.clean_error( );  

    Data.Solve( Error, Msg );

//    try{
//      Data.CheckResidualConvergence( Error , Msg );
//    } catch ( MAP_ERROR_CODE code ) { 
//      std::string str = "";                                   
//      str += boost::lexical_cast < std::string > ( code );    
//      if ( str.size( ) == 1 ) str = " " + str;                
//      str += "] : " + MAP_ERROR_CODE_to_string.at( code );    
//
//      Error.set_error_key( MAP_ERROR );                       
//      Msg.RecordErrorToErrorList( str );                                  
//    }

    // this sets the initial conditions for the Fix and Vessel nodes
    // This adds sum_FX, sum_FY and sum_FZ into the node
    for( int i=0 ; i<Data.getSizeOfNode() ; i++) {
      Data.initializeCableNode( i );
    };//END for
  }
  else {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_54 );
    str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_54 );

    Msg.RecordErrorToErrorList( str );
    Error.set_error_key( MAP_ERROR );  
  }//END if
};


/**
 * ====================================================================================================
 * MSQS_CalcOutput
 * 
 * The following flags are evaluated:
 *     - PLOT_flag        
 *     - X_POS_flag       
 *     - Y_POS_flag       
 *     - Z_POS_flag       
 *     - X_FORCE_flag     
 *     - Y_FORCE_flag     
 *     - Z_FORCE_flag     
 *     - LINE_TENSION_flag
 *     - OMIT_CONTACT_flag
 *     - LAY_LENGTH_flag  
 *
 * ====================================================================================================
 */
void MSQS_CalcOutput( float T                          ,
		      MAP_InputType_class           &In      ,
		      MAP_ParameterType_class       &Param   ,
		      void* NULL_ContinuousState       ,
		      void* NULL_DiscreteState         ,
		      MAP_ConstraintStateType_class &Constrnt,
		      MAP_OtherStateType_class      &Data    ,
		      MAP_OutputType_class          &Out     ,
		      MAP_ErrStat             &Error   ,
		      MAP_Message             &Msg ){

  Msg.messageClean();   // clean the Msg parameter
    
  if ( Data.isNumericsUninitialized()==false && Error.error_status()==MAP_SAFE ){

    Error.clean_error( ); // set MAP error status (presuming here there is no error
                          // just by simply extering this function from the glue code)    

    if ( Data.isOutputFileOpen() ){
      for(int i=0 ; i<Data.getSizeOfElement() ; i++ ){
        Data.getOutputStreamValue( i , T , Msg );
      };//END for
      
      Data.writeToOutputFile();	
      Data.cleanFileStreamBuffer();
    }
    else {
      // open file we are writting the output to
      Data.openOutputFile();
      
      // Get current date/time, format is YYYY-MM-DD.HH:mm:ss
      time_t     now = time(0);
      struct tm  tstruct;
      char       buf[80];
      tstruct = *localtime(&now);
      
      strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
      
      Data.writeToOutputString("MAP Output File.\n");
      Data.writeToOutputString("Outputs were generated using MAP v1.0 on <");
      Data.writeToOutputString(buf);
      Data.writeToOutputString(">\n\n");

      // print the line telling user what each column represents
      for(int i=0 ; i<Data.getSizeOfElement() ; i++ ){
        Data.getOutputStreamHeader( i , Msg );
      };//END for

      // print the units for each column
      for(int i=0 ; i<Data.getSizeOfElement() ; i++ ){
        Data.getOutputStreamUnits( i , Msg );
      };//END for

      Data.writeToOutputFile();	
      Data.cleanFileStreamBuffer();

      for(int i=0 ; i<Data.getSizeOfElement() ; i++ ){
        Data.getOutputStreamValue( i , T , Msg );
      };//END for

      Data.writeToOutputFile();	
      Data.cleanFileStreamBuffer();
    };//END if
  }
  else {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_68 );
    str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_68 );

    Msg.RecordErrorToErrorList( str );
    Error.set_error_key( MAP_ERROR );  
  }//END if

  Data.plot(Error,Msg);
};


/**
 * ====================================================================================================
 * MSQS_End( )
 * ====================================================================================================
 */
void MSQS_End( MAP_InputType_class           &In      ,
               MAP_ParameterType_class       &Param   ,
               void* NULL_ContinuousState       ,
               void* NULL_DiscreteState         ,
               MAP_ConstraintStateType_class &Constrnt,
               MAP_OtherStateType_class      &Data    ,
               MAP_OutputType_class          &Out     ,
               MAP_ErrStat             &Error   ,
               MAP_Message             &Msg ){

  Error.clean_error( ); // set MAP error status (presuming here there is no error
                        // just by simply extering this function from the glue code)    
  Msg.messageClean();

  if ( Data.isNumericsUninitialized()==false ){
    // free all PETSc data
    // @todo : a seg-fault is created when:
    //     1) Init.setGravity is commented in the python
    //     2) the node type is misspelled
    Data.cleanNumericSolver( Error, Msg );
    Data.closeOutputFile( );
  }
  else{	
    std::string str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_65 );
    str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_65 );

    Msg.RecordErrorToErrorList( str );
    Error.set_error_key( MAP_ERROR );  
  };//END if
};
