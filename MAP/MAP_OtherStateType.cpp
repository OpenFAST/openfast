/**
 * ====================================================================================================
 *                              MAP_OtherStateType_class.cpp
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


#include "MAP_OtherStateType.h"


/**
 * ============================================================================
 * addCableLibrary
 *
 * Adds a cable type (such as steel, nylon, vectran) to the cable library
 *
 * @input : T[0]         -- name of variable
 * @input : T[1]         -- diameter
 * @input : T[2]         -- density of material
 * @input : T[3]         -- Young's modulus
 * @input : T[4]         -- friction coef. between seabed and cable 
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void MAP_OtherStateType_class::addCableLibrary( const std::vector<std::string> &T     ,
                                          const int                      index  ,
                                          MAP_ErrStat                    &Error ,
                                          MAP_Message                    &Msg ) {
    
  CableLibrary_ptr prop_ptr( new CableLibrary );      // create new memory on the heap (dynamically allocated)

  prop_ptr->label     = T[0];                           // First give the cable a name, such as steel or nylon

  /**
   * =======  Set cable diameter properties  ======     <------------------------------------------------------+
   */                                                                                              //          |
  //          |
  prop_ptr->Diam.name = "Diam";                                                                    //          |
  //          |
  if ( boost::starts_with( T[1] , "#" ) ) { // if the first elemt in diameter starts with '#'      //          |
    Error.set_error_key( MAP_WARNING );   // This is a warning printed to the message string     //          |
    Msg.RecordWarningToWarningList("Ignoring '#' chatacter preceeding 'Diam' in CableLibrary parameter.");     //          |
    //          |
    // Now cycle throw the string and write it to the This parameter,                            //          |
    // but ignore the first character (which will be a '#' due to                                //          |
    // if-statement condition above                                                              //          |
    std::string This = "";                                                                       //          |
    for ( unsigned int i=1 ; i<T[1].size() ; i++) {                                              //          |
      This.push_back( T[1][i] );                                                               //          |
    }                                                                                            //          |
    //          |
    try { // make sure the variable can be converted to a double                                 //          |
      prop_ptr->Diam.value = boost::lexical_cast<double>( This );                              //          |
      prop_ptr->Diam.index = index;
    } catch ( boost::bad_lexical_cast& ) {                                                       //          |
      std::string str = " ";                                                                   //          |
      str += boost::lexical_cast < std::string > ( MAP_ERROR_6 );                              //          |
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_6 );                              //          |
      //          |
      Error.set_error_key( MAP_ERROR );                                                        //          |
      //          |
      // Here is the exception message                                                         //          |
      Msg.RecordErrorToErrorList( str );                                                                   //          |
      //          |
      // create context to write error in summary.map                                          //          |
      str.erase(2,1);                                                                          //          |
      std::ostringstream S;                                                                    //          |
      S.str("");S.clear();                                                                     //          |
      S << ">>>> " +str +"    :    MAP error in file ";                                        //          |
      S << __FILE__;                                                                           //          |
      S << " on line ";                                                                        //          |
      S << __LINE__;                                                                           //          |
      S << " occurred.";                                                                       //          |
      Msg.WriteErrorToOutputFile( S.str() );                                                              //          |
    };//END try                                                                                  //          |
  }                                                                                                //          |
  else {                                                                                           //          |
    try { // make sure the variable can be converted to a double                                 //          |
      prop_ptr->Diam.value = boost::lexical_cast<double>( T[1] );                              //          |
      prop_ptr->Diam.index = index;
    } catch ( const boost::bad_lexical_cast& ) {                                                 //          |
      std::string str = " ";                                                                   //          |
      str += boost::lexical_cast < std::string > ( MAP_ERROR_6 );                              //          |
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_6 );                              //          |
      //          |
      Error.set_error_key( MAP_ERROR );                                                        //          |
      //          |
      // Here is the exception message                                                         //          |
      Msg.RecordErrorToErrorList( str );                                                                   //          |
      //          |
      // create context to write error in summary.map                                          //          |
      str.erase(2,1);                                                                          //          |
      std::ostringstream S;                                                                    //          |
      S.str("");S.clear();                                                                     //          |
      S << ">>>> " +str +"    :    MAP error in file ";                                        //          |
      S << __FILE__;                                                                           //          |
      S << " on line ";                                                                        //          |
      S << __LINE__;                                                                           //          |
      S << " occurred.";                                                                       //          |
      Msg.WriteErrorToOutputFile( S.str() );                                                              //          |
    };//END try                                                                                  //          |
  };                                                                                               // ---------+
  //============== <END> =======================================================================================
    
  /**
   * =======  Set cable density properties  ======     <-------------------------------------------------------+
   */                                                                                              //          |
  //          |
  prop_ptr->MassDenInAir.name = "MassDenInAir";                 // set density properties                  //          |
  //          |
  if ( boost::starts_with( T[2] , "#" ) ) {  // if the first character in MassDenInAir starts with '#' //          |
    //          |
    Error.set_error_key( MAP_WARNING ); // This is a warning printed to the message string       //          |
    Msg.RecordWarningToWarningList("Ignoring '#' chatacter preceeding 'MassDenInAir' in CableLibrary parameter."); //          |
                                                                                                     //          |
    // Now cycle throw the string and write it to the This parameter,                            //          |
    // but ignore the first character (which will be a '#' due to                                //          |
    // if-statement condition above                                                              //          |
    std::string This = "";                                                                       //          |
    for ( unsigned int i=1 ; i<T[2].size() ; i++) {                                              //          |
      This.push_back( T[2][i] );                                                               //          |
    };// END for                                                                                 //          |
    //          |
    try { // make sure the variable can be converted to a double                                 //          |
      prop_ptr->MassDenInAir.value = boost::lexical_cast<double>( This );                      //          |
      prop_ptr->MassDenInAir.index = index;
    } catch ( boost::bad_lexical_cast& ) {                                                       //          |
      std::string str = "";                                                                    //          |
      str += boost::lexical_cast < std::string > ( MAP_ERROR_17 );                             //          |
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_17 );                             //          |
      //          |
      Error.set_error_key( MAP_ERROR );                                                        //          |
      //          |
      // Here is the exception message                                                         //          |
      Msg.RecordErrorToErrorList( str );                                                                   //          |
      //          |
      // create context to write error in summary.map                                          //          |
      str.erase(2,1);                                                                          //          |
      std::ostringstream S;                                                                    //          |
      S.str("");S.clear();                                                                     //          |
      S << ">>>> " +str +"    :    MAP error in file ";                                        //          |
      S << __FILE__;                                                                           //          |
      S << " on line ";                                                                        //          |
      S << __LINE__;                                                                           //          |
      S << " occurred.";                                                                       //          |
      Msg.WriteErrorToOutputFile( S.str() );                                                              //          |
    };// END try                                                                                 //          |
  }                                                                                                //          |
  else {                                                                                           //          |
    try { // make sure the variable can be converted to a double                                 //          |
      prop_ptr->MassDenInAir.value = boost::lexical_cast<double>( T[2] );                      //          |
      prop_ptr->MassDenInAir.index = index;                                                    //          |
    } catch ( boost::bad_lexical_cast& ) {                                                       //          |
      std::string str = "";                                                                    //          |
      str += boost::lexical_cast < std::string > ( MAP_ERROR_17 );                             //          |
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_17 );                             //          |
      //          |
      Error.set_error_key( MAP_ERROR );                                                        //          |
      //          |
      // Here is the exception message                                                         //          |
      Msg.RecordErrorToErrorList( str );                                                                   //          |
      //          |
      // create context to write error in summary.map                                          //          |
      str.erase(2,1);                                                                          //          |
      std::ostringstream S;                                                                    //          |
      S.str("");S.clear();                                                                     //          |
      S << ">>>> " +str +"    :    MAP error in file ";                                        //          |
      S << __FILE__;                                                                           //          |
      S << " on line ";                                                                        //          |
      S << __LINE__;                                                                           //          |
      S << " occurred.";                                                                       //          |
      Msg.WriteErrorToOutputFile( S.str() );                                                              //          |
    };// END try                                                                                 //          |
  };                                                                                               // ---------+
  //============== <END> =======================================================================================


  /**
   * =======  Set cable Young's modulus properties  ======     <-----------------------------------------------+
   */                                                                                              //          |
  //          |
  prop_ptr->EA.name = "EA";                                                                        //          |
  //          |
  if ( boost::starts_with( T[3] , "#" ) ) { // if the first character in 'E' starts with '#'       //          |
    //          |
    Error.set_error_key( MAP_WARNING ); // This is a warning printed to the message string       //          |
    Msg.RecordWarningToWarningList("Ignoring '#' chatacter preceeding 'EA' in CableLibrary parameter.");       //          |
    //          |
    // Now cycle throw the string and write it to the This parameter,                            //          |
    // but ignore the first character (which will be a '#' due to                                //          |
    // if-statement condition above                                                              //          |
    std::string This = "";                                                                       //          |
    for ( unsigned int i=1 ; i<T[3].size() ; i++) {                                              //          |
      This.push_back( T[3][i] );                                                               //          |
    };// END for                                                                                 //          |
    //          |
    try { // make sure the variable can be converted to a double                                 //          |
      prop_ptr->EA.value = boost::lexical_cast<double>( This );                                //          |
      prop_ptr->EA.index = index;
    } catch ( boost::bad_lexical_cast& ) {                                                       //          |
      std::string str = "";                                                                    //          |
      str += boost::lexical_cast < std::string > ( MAP_ERROR_18 );                             //          |
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_18 );                             //          |
      //          |
      Error.set_error_key( MAP_ERROR );                                                        //          |
      //          |
      // Here is the exception message                                                         //          |
      Msg.RecordErrorToErrorList( str );                                                                   //          |
      //          |
      // create context to write error in summary.map                                          //          |
      str.erase(2,1);                                                                          //          |
      std::ostringstream S;                                                                    //          |
      S.str("");S.clear();                                                                     //          |
      S << ">>>> " +str +"    :    MAP error in file ";                                        //          |
      S << __FILE__;                                                                           //          |
      S << " on line ";                                                                        //          |
      S << __LINE__;                                                                           //          |
      S << " occurred.";                                                                       //          |
      Msg.WriteErrorToOutputFile( S.str() );                                                              //          |
    };// END try                                                                                 //          |
  }                                                                                                //          |
  else {                                                                                           //          |
    try { // make sure the variable can be converted to a double                                 //          |
      prop_ptr->EA.value = boost::lexical_cast<double>( T[3] );                                //          |
      prop_ptr->EA.index = index;
    } catch ( boost::bad_lexical_cast& ) {                                                       //          |
      std::string str = "";                                                                    //          |
      str += boost::lexical_cast < std::string > ( MAP_ERROR_18 );                             //          |
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_18 );                             //          |
      //          |
      Error.set_error_key( MAP_ERROR );                                                        //          |
      //          |
      // Here is the exception message                                                         //          |
      Msg.RecordErrorToErrorList( str );                                                                   //          |
      //          |
      // create context to write error in summary.map                                          //          |
      str.erase(2,1);                                                                          //          |
      std::ostringstream S;                                                                    //          |
      S.str("");S.clear();                                                                     //          |
      S << ">>>> " +str +"    :    MAP error in file ";                                        //          |
      S << __FILE__;                                                                           //          |
      S << " on line ";                                                                        //          |
      S << __LINE__;                                                                           //          |
      S << " occurred.";                                                                       //          |
      Msg.WriteErrorToOutputFile( S.str() );                                                              //          |
    };// END try                                                                                 //          |
  };                                                                                               // ---------+
  //============== <END> =======================================================================================


  /**
   * =======  Set cable/sea bed friction coefficients  ======     <--------------------------------------------+
   */                                                                                              //          |
  //          |
  prop_ptr->CB.name = "CB";                                                                        //          |
  //          |
  if ( boost::starts_with( T[4] , "#" ) ) { // if the first character in 'CB starts with '#'       //          |
    //          |
    Error.set_error_key( MAP_WARNING ); // This is a warning printed to the message string       //          |
    Msg.RecordWarningToWarningList("Ignoring '#' chatacter preceeding 'CB' in CableLibrary parameter.");       //          |
    //          |
    // Now cycle throw the string and write it to the This parameter,                            //          |
    // but ignore the first character (which will be a '#' due to                                //          |
    // if-statement condition above                                                              //          |
    std::string This = "";                                                                       //          |
    for ( unsigned int i=1 ; i<T[4].size() ; i++) {                                              //          |
      This.push_back( T[4][i] );                                                               //          |
    };// END for                                                                                 //          |
    //          |
    try { // make sure the variable can be converted to a double                                 //          |
      prop_ptr->CB.value = boost::lexical_cast<double>( This );                                //          |
      prop_ptr->CB.index = index;
    } catch ( boost::bad_lexical_cast& ) {                                                       //          |
      std::string str = "";                                                                    //          |
      str += boost::lexical_cast < std::string > ( MAP_ERROR_19 );                             //          |
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_19 );                             //          |
      //          |
      Error.set_error_key( MAP_ERROR );                                                        //          |
      //          |
      // Here is the exception message                                                         //          |
      Msg.RecordErrorToErrorList( str );                                                                   //          |
      //          |
      // create context to write error in summary.map                                          //          |
      str.erase(2,1);                                                                          //          |
      std::ostringstream S;                                                                    //          |
      S.str("");S.clear();                                                                     //          |
      S << ">>>> " +str +"    :    MAP error in file ";                                        //          |
      S << __FILE__;                                                                           //          |
      S << " on line ";                                                                        //          |
      S << __LINE__;                                                                           //          |
      S << " occurred.";                                                                       //          |
      Msg.WriteErrorToOutputFile( S.str() );                                                              //          |
    };// END try                                                                                 //          |
  }                                                                                                //          |
  else {                                                                                           //          |
    try { // make sure the variable can be converted to a double                                 //          |
      prop_ptr->CB.value = boost::lexical_cast<double>( T[4] );                                //          |
      prop_ptr->CB.index = index;
    } catch ( boost::bad_lexical_cast& ) {                                                       //          |
      std::string str = "";                                                                    //          |
      str += boost::lexical_cast < std::string > ( MAP_ERROR_19 );                             //          |
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_19 );                             //          |
      //          |
      Error.set_error_key( MAP_ERROR );                                                        //          |
      //          |
      // Here is the exception message                                                         //          |
      Msg.RecordErrorToErrorList( str );                                                                   //          |
      //          |
      // create context to write error in summary.map                                          //          |
      str.erase(2,1);                                                                          //          |
      std::ostringstream S;                                                                    //          |
      S.str("");S.clear();                                                                     //          |
      S << ">>>> " +str +"    :    MAP error in file ";                                        //          |
      S << __FILE__;                                                                           //          |
      S << " on line ";                                                                        //          |
      S << __LINE__;                                                                           //          |
      S << " occurred.";                                                                       //          |
      Msg.WriteErrorToOutputFile( S.str() );                                                              //          |
    };// END try                                                                                 //          |
  };                                                                                               // ---------+
  //============== <END> =======================================================================================

  // push the new object we created into our vector container CableLibrary
  property.push_back( prop_ptr );
};


/**
 * ============================================================================
 * writeCableLibraryData
 *
 * write all data from the CableLibrary class to the MAP_Message parameter
 *
 * @input : Msg          -- Error message status
 * ============================================================================
 */
void MAP_OtherStateType_class::writeCableLibraryData( MAP_Message &Msg) {

    
  for ( unsigned int i=0 ; i<property.size() ; i++){

    std::ostringstream S;
    S << std::fixed << std::setprecision(3);
	
    /**
     * Line 1 : 
     *
     * write the cable type to the Msg parameter
     */
    Msg.WriteDataToOutputFile("    Cable Type          : ");
    Msg.WriteDataToOutputFile( property[i]->label );
    Msg.WriteDataToOutputFile("\n");
	    
    /**
     * Line 2 : 
     *
     * write the cable diameter to the Msg parameter
     */
    S << property[i]->Diam.value;
    Msg.WriteDataToOutputFile("    ");
    Msg.WriteDataToOutputFile( property[i]->Diam.name );
    // Msg.WriteDataToOutputFile(":     ");
    Msg.WriteDataToOutputFile("         [m]    : ");
    Msg.WriteDataToOutputFile( S.str() );
    Msg.WriteDataToOutputFile("\n");    
    S.str("");S.clear();

    /**
     * Line 3 : 
     *
     * write the cable density to the Msg parameter
     */
    S << property[i]->MassDenInAir.value;
    Msg.WriteDataToOutputFile("    ");
    Msg.WriteDataToOutputFile( property[i]->MassDenInAir.name );
    Msg.WriteDataToOutputFile(" [kg/m] : ");
    Msg.WriteDataToOutputFile( S.str() );
    Msg.WriteDataToOutputFile("\n");    
    S.str("");S.clear();

    /**
     * Line 5 : 
     *
     * write the cable Young's modulus to the Msg parameter
     */
    S << property[i]->EA.value;
    Msg.WriteDataToOutputFile("    ");
    Msg.WriteDataToOutputFile( property[i]->EA.name );
    Msg.WriteDataToOutputFile("           [N]    : ");
    Msg.WriteDataToOutputFile( S.str() );
    Msg.WriteDataToOutputFile("\n");    
    S.str("");S.clear();

    /**
     * Line 6 : 
     *
     * write the cable friction coefficient to the Msg parameter
     */
    S << property[i]->CB.value;
    Msg.WriteDataToOutputFile("    ");
    Msg.WriteDataToOutputFile( property[i]->CB.name );
    Msg.WriteDataToOutputFile("                  : ");
    Msg.WriteDataToOutputFile( S.str() );
    Msg.WriteDataToOutputFile("\n\n");    
    S.str("");S.clear();
  }//END for
};


/**
 * ============================================================================
 * writeEnvironmentData
 *
 * prints the settings for the sea density, water depth and gravity constrant
 * 
 * @input : Msg          -- Error message status
 * ============================================================================
 */
void MAP_OtherStateType_class::writeEnvironmentData  ( MAP_Message &Msg ){
  std::ostringstream S;
  S << std::fixed << std::setprecision(3);

  // Line 1: write sea density to the Msg parameter
  S << rho_sea;
  Msg.WriteDataToOutputFile("    Sea Density [kg/m^3]: ");
  Msg.WriteDataToOutputFile( S.str() );
  Msg.WriteDataToOutputFile("\n");    
  S.str("");S.clear();
	
  // Line 2: write gravity  to the Msg parameter
  S << gravity;
  Msg.WriteDataToOutputFile("    Gravity     [m/s^2] : ");
  Msg.WriteDataToOutputFile( S.str() );
  Msg.WriteDataToOutputFile("\n");    
  S.str("");S.clear();
	
  // Line 3: write depth to the Msg parameter
  S << depth;
  Msg.WriteDataToOutputFile("    Depth       [m]     : ");
  Msg.WriteDataToOutputFile( S.str() );
  Msg.WriteDataToOutputFile("\n\n");    
  S.str("");S.clear();
};


/**
 * ============================================================================
 * writeNodeData
 *
 * @input : Msg          -- Error message status
 * ============================================================================
 */
void MAP_OtherStateType_class::writeNodeData( MAP_Message &Msg ){
  std::vector < std::vector<std::string> > cell;
  std::string write = "";
  std::ostringstream S;
  S << std::fixed << std::setprecision(3);

  for(int row=0 ; row<10 ; row++) {
    // This adds a row
    cell.push_back( std::vector<std::string>() );
  }// END for	
    
  // This creates four columns of node data per line.
  unsigned int FOUR = 4;

  for (unsigned int j=0 ; j<node.size() ; j=j+FOUR) {
    int num = 0;

    if( j+FOUR > node.size() ){
      num = node.size() - j;
    }
    else{
      num = FOUR;
    };//END if
	
    for ( unsigned int col=j ; col<j+num ; col++) {
	    
      // This creates the header for the Node data portion of the output message
      write.clear();
	    
      S << col+1;
      // write += "   Node ";
      write += "Node ";
      write += S.str();
      write += " Data:";
      S.str("");S.clear();
	    
      /**
       * To be sure the length of 'write' is not greater than
       * length we are allow to across, this assert statement 
       * is raised.
       */
      assert( write.size()<_TEXT_COLOR::STR_LEN );
	    
      /**
       * Fill spaces in the 'write' string with white spaces to 
       * occupy area in the output string (to make the output
       * print in alignment)
       */
      do {
        write += " ";
      } while( write.size()<_TEXT_COLOR::STR_LEN);

      // Add column to all rows
      cell[0].push_back( write );
      write.clear();

      switch ( node[col]->getNodeType() ) {
      case No_Definition:
        // write += "   No Definition ";
        write += "No Definition ";
        break;

      case Fix : 
        // write += "   Fix";
        write += "Fix";
        break;
	    
      case Vessel :
        // write += "   Vessel";
        write += "Vessel";
        break;
	    
      case Connect : 
        // write += "   Connect" ;
        write += "Connect";
        break;
      };// END switch

      /**
       * To be sure the length of 'write' is not greater than
       * length we are allow to across, this assert statement 
       * is raised. This only because visiable if the condition
       * is violated
       */
      assert( write.size()<_TEXT_COLOR::STR_LEN );

      /**
       * Fill spaces in the 'write' string with white spaces to 
       * occupy area in the output string (to make the output
       * print in alignment)
       */
      do{
        write += " ";
      } while( write.size()<_TEXT_COLOR::STR_LEN );
	
      // Add column to all rows
      cell[1].push_back( write );
      write.clear();
	
      /**
       * this block pushes strings to each cell in the output string
       * we are writting to
       */
      cell[2].push_back( VarType::writeNodeData( node[col]->X ) );
      cell[3].push_back( VarType::writeNodeData( node[col]->Y ) );
      cell[4].push_back( VarType::writeNodeData( node[col]->Z ) );
      cell[5].push_back( VarType::writeNodeData( node[col]->M ) );
      cell[6].push_back( VarType::writeNodeData( node[col]->B ) );
      cell[7].push_back( VarType::writeNodeData( node[col]->FX ) );
      cell[8].push_back( VarType::writeNodeData( node[col]->FY ) );
      cell[9].push_back( VarType::writeNodeData( node[col]->FZ ) );
    }// END for

    // Write node number to message string
    write += "            |  " ;
    for( unsigned int i=j ; i<j+num ; i++) write += cell[0][i]; 
    write += "\n";
    
    write += "            |  ---------------------------------------------------------------------------------------\n" ;

    // Write node type to message string
    write += "Node Type:  |  " ;
    for( unsigned int i=j ; i<j+num ; i++) write += cell[1][i]; 
    write += "\n";
  
    // Write X data to message string
    write += "X  [m]   :  |  " ;
    for( unsigned int i=j ; i<j+num ; i++) write += cell[2][i]; 
    write += "\n";
  
    // Write Y data to message string
    write += "Y  [m]   :  |  " ;
    for( unsigned int i=j ; i<j+num ; i++) write += cell[3][i]; 
    write += "\n";
  
    // Write Z data to message string
    write += "Z  [m]   :  |  " ;
    for( unsigned int i=j ; i<j+num ; i++) write += cell[4][i]; 
    write += "\n";
  
    // Write M data to message string
    write += "M  [kg]  :  |  " ;
    for( unsigned int i=j ; i<j+num ; i++) write += cell[5][i]; 
    write += "\n";
  
    // Write B data to message string
    write += "B  [m^3] :  |  " ;
    for( unsigned int i=j ; i<j+num ; i++) write += cell[6][i]; 
    write += "\n";

    // Write FX data to message string
    write += "FX [N]   :  |  " ;
    for( unsigned int i=j ; i<j+num ; i++) write += cell[7][i];
    write += "\n";

    // Write FY data to message string    
    write += "FY [N]   :  |  " ;
    for( unsigned int i=j ; i<j+num ; i++) write += cell[8][i]; 
    write += "\n";
  
    // Write FZ data to message string
    write += "FZ [N]   :  |  " ;
    for( unsigned int i=j ; i<j+num ; i++) write += cell[9][i]; 
    write += "\n\n\n";

    // Finalize output by writting the string to the Msg object
    Msg.WriteDataToOutputFile( write );
  }
};


/**
 * ============================================================================
 * addNode  >>  Overloaded Member Function
 *
 * @input : T[0]         --  name of variable
 * @input : T[1]         --  X
 * @input : T[2]         --  Y
 * @input : T[3]         --  Z
 * @input : T[4]         --  M
 * @input : T[5]         --  B
 * @input : T[6]         --  FX
 * @input : T[8]         --  FY
 * @input : T[9]         --  FZ
 * @input : index        -- 
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void MAP_OtherStateType_class::addNode(const std::vector<std::string> &T     ,
                                 const int                      index  ,
                                 MAP_ErrStat                    &Error ,
                                 MAP_Message                    &Msg) {
    
  // throw an exception if the node type in the MAP input file
  // does not match one of the specified types
  // create new memory on the heap (dynamically allocated)

  NodeType node_type;

  if     ( boost::iequals( T[1] , "FIX"     ) ) node_type = Fix;      // Fix node = 1 enum type
  else if( boost::iequals( T[1] , "CONNECT" ) ) node_type = Connect;  // Connect node = 2 enum type 
  else if( boost::iequals( T[1] , "VESSEL"  ) ) node_type = Vessel;   // Vessel node = 3 enum type
  else throw MAP_ERROR_15;                                            // exception thrown otherwise

  Node_ptr node_ptr( new Node( node_type , this->gravity , this->rho_sea ) );

  // set node parameters
  node_ptr->setVarType( T[2] , "X"  , index , &Node::X , Error, Msg );
  node_ptr->setVarType( T[3] , "Y"  , index , &Node::Y , Error, Msg );

  // the case for "Z" is handled a few lines below
  node_ptr->setVarType( T[5] , "M"  , index , &Node::M , Error, Msg );
  node_ptr->setVarType( T[6] , "B"  , index , &Node::B , Error, Msg );
  node_ptr->setVarType( T[7] , "FX" , index , &Node::FX, Error, Msg );
  node_ptr->setVarType( T[8] , "FY" , index , &Node::FY, Error, Msg );
  node_ptr->setVarType( T[9] , "FZ" , index , &Node::FZ, Error, Msg );

  // This if-statement evaluates whether the Z variable in the MAP input file reads 
  // as 'depth'. If it does, then we use the depth marked in the calling program. 
  if( boost::iequals( T[4] , "DEPTH" ) ){
    // set the Z node displacement to the default water depth
    node_ptr->setVarType( boost::lexical_cast<std::string>( this->getDepth() ) , "Z" , index , &Node::Z , Error , Msg );
  }
  else{
    node_ptr->setVarType( T[4] , "Z"  , index , &Node::Z , Error, Msg );
  };// END if

    // Now psuh the node we just created to the stack
  node.push_back( node_ptr );	

  // Check if Z is smaller than the water depth. This cannot happen. Warn the users they
  // used a potentially wrong value in the MAP input file.
  if ( node_ptr->getPtrValue( &Node::Z ) < this->getDepth() ){
    node_ptr->setPtrValue( &Node::Z , this->getDepth() );

    // issue error
    throw MAP_ERROR_16;
  };// END if
};


/**
 * ============================================================================
 * addElement  
 * ============================================================================
 */
void MAP_OtherStateType_class::addElement( const std::vector<std::string> &T     ,
                                     const int                      index  ,
                                     MAP_ErrStat                    &Error ,
                                     MAP_Message                    &Msg ) {

  // create new memory on the heap (dynamically allocated)
  Element_ptr element_ptr( new Element );
  std::string error_output = "";

  element_ptr->setLineProperty( T[1] , property, Error, Msg);
  
  // set fairlead node 
  //  
  // Sets the properties for the bottom node.  If the input file is not formatted 
  // correctly, an exception is thrown. -9999 is a default flag notifying MAP the  
  // parameter is uninitialized. In this case, if the parameter is uninitialized, 
  // the program will fail. This is considered a critical failure, and needs 
  // to be corrected before any calculations can be performed.
  unsigned int top    = 9999;
  unsigned int bottom = 9999;
  try { // make sure the variable can be converted to a double
    top = boost::lexical_cast<int>( T[4] );

    // make sure 'top' is not out of index
    if ( top > node.size( ) ) throw MAP_ERROR_21;

    // 'node' is a vector of nodes. 'top' is  the index we are setting
    element_ptr->setFairlead( top , node );

    // convert string to integer
    bottom = boost::lexical_cast<int>( T[3] );

    // make sure 'top' is not out of index
    if ( bottom > node.size( ) ) throw MAP_ERROR_22;

    // 'node' is a vector of nodes. 'bottom' is
    // the index we are setting
    element_ptr->setAnchor  ( bottom , node );
    
    // set element unstretched length
    // Here we set the Lu variable in class Element
    element_ptr->setLu( T[2], index, Error, Msg );

    // Finally, set H and V as iterated value
    element_ptr->setH_and_V_flags(); 
    element_ptr->setH( index , Error , Msg ); // this calls VarType::setGenericVarType
    element_ptr->setV( index , Error , Msg ); // this calls VarType::setGenericVarType

    element_ptr->setHX(); // purely for writing Hx data to the output file
    element_ptr->setHY(); // purely for writing Hy data to the output file

    int cnt = 0;
    if( element_ptr->getLuFlag()==false ) cnt++;
    if( element_ptr->getHFlag() ==false ) cnt++;
    if( element_ptr->getVFlag() ==false ) cnt++;

    // if too many variables are iterated, issue a warning
    if (cnt >=3 ) throw MAP_ERROR_24;

    // if too little variables are iterated, issue a warning
    if (cnt <=1 ) throw MAP_ERROR_25;

    for (unsigned int i=5 ; i<T.size()-1 ; i++ ){
      this->setElementOptions( element_ptr , T[i] , Error, Msg );
    };// END for

    // This last step creates the node. This is a boost::shared_ptr
    element.push_back( element_ptr );
  } catch ( boost::bad_lexical_cast const& ) {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_23 );
    str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_23 );

    Error.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordErrorToErrorList( str );

    // create context to write error in summary.map     
    str.erase(2,1);                                     
    std::ostringstream S;                               
    S.str("");S.clear();                                
    S << ">>>> " +str +"    :    MAP error in file ";   
    S << __FILE__;                                      
    S << " on line ";                                   
    S << __LINE__;                                      
    S << " occurred.";                                  
    Msg.WriteErrorToOutputFile( S.str() );                         
  } catch ( MAP_ERROR_CODE &code ) {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( code );
    std::ostringstream S;
    S << index+1;
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str( ) + ".";

    Error.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordErrorToErrorList( str );

    // create context to write error in summary.map        
    str.erase(2,1);                                        
    S.str("");S.clear();                                   
    S << ">>>> " +str +"    :    MAP error in file ";      
    S << __FILE__;                                         
    S << " on line ";                                      
    S << __LINE__;                                         
    S << " occurred.";                                     
    Msg.WriteErrorToOutputFile( S.str() );                            
  };//END try
};


/**
 * ============================================================================
 * associateElementVarTypeToNWTCType
 *
 * @input : i            -- 
 * @input : P            -- 
 * @input : C            --   
 *
 * @todo : V and H should be local states instead of assigned to parameters
 *         (since they change in time)
 * ============================================================================
 */
void MAP_OtherStateType_class::associateElementVarTypeToNWTCType( const int               index ,
                                                            MAP_ParameterType_class       &P    , 
                                                            MAP_ConstraintStateType_class &C    ,
                                                            MAP_ErrStat             &Err  ,
                                                            MAP_Message             &Msg) {
  // Assign element H variable
  if ( this->elementVarTypeBool( index , &Element::H )==true ) {
    this->setElementReferenceToMAPType( index , P , &Element::H );  
  }   
  else {                                                                    
    this->setElementReferenceToMAPType( index , C , &Element::H );  
  };                                                                        
                                                                               	
  // Assign element V variable
  if ( this->elementVarTypeBool( index , &Element::V )==true ) {              
    this->setElementReferenceToMAPType( index , P , &Element::V );     
  }
  else {                                                                    
    this->setElementReferenceToMAPType( index , C , &Element::V );  
  };                                                                        
                                                                              
  // Assign element Lu variable
  if ( this->elementVarTypeBool( index , &Element::Lu )==true ){ 
    this->setElementReferenceToMAPType( index , P , &Element::Lu );    
  }
  else {                                                                    
    this->setElementReferenceToMAPType( index , C , &Element::Lu ); 
  };
};


/**
 * ============================================================================
 * checkElementVarTypeReferences
 *
 * As a last step, we check the reference_counter for all VarType (in Node 
 * Element and CableLibrary) to ensure it is pointing to one MAP_ParameterType_class,               
 * MAP_ConstraintStateType_class or MAP_InputType_class.
 *
 * @input : i            -- 
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void MAP_OtherStateType_class::checkElementVarTypeReferences( const int   index  , 
							MAP_ErrStat &Error , 
							MAP_Message &Msg ){
  if ( this->checkElementReference( index , &Element::Lu ) != 1 ){    
    std::string str = "";
    std::ostringstream S;
    S << index+1;
    str += boost::lexical_cast < std::string > ( MAP_ERROR_42 );
    str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_42 ) + S.str() + ".";

    Error.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordErrorToErrorList( str );

  };// END if
};


/**
 * ============================================================================
 * checkNodeVarTypeReferences 
 *
 * As a last step, we check the reference_counter for all VarType (in Node 
 * Element and CableLibrary) to ensure it is pointing to one MAP_ParameterType_class,               
 * MAP_ConstraintStateType_class or MAP_InputType_class.
 *
 * @input : i            -- 
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void MAP_OtherStateType_class::checkNodeVarTypeReferences( const int index    , 
						     MAP_ErrStat &Error , 
						     MAP_Message &Msg ){

  try {
    if ( this->checkNodeReference( index , &Node::X  ) != 1 ){
      throw MAP_ERROR_43;
    };// END if

    if ( this->checkNodeReference( index , &Node::Y  ) != 1 ){
      throw MAP_ERROR_44;
    };// END if

    if ( this->checkNodeReference( index , &Node::Z  ) != 1 ){    
      throw MAP_ERROR_45;
    };// END if

    if ( this->checkNodeReference( index , &Node::M  ) != 1 ){    
      throw MAP_ERROR_46;
    };// END if

    if ( this->checkNodeReference( index , &Node::B  ) != 1 ){    
      throw MAP_ERROR_47;
    };// END if

    if ( this->checkNodeReference( index , &Node::FX ) != 1 ){   
      throw MAP_ERROR_48;
    };// END if

    if ( this->checkNodeReference( index , &Node::FY ) != 1 ){    
      throw MAP_ERROR_49;
    };// END if

    if ( this->checkNodeReference( index , &Node::FZ ) != 1 ){    
      throw MAP_ERROR_50;
    };// END if
  } catch( MAP_ERROR_CODE &code ) {
    std::string str = "";
    std::ostringstream S;
    S << index+1;
    str += boost::lexical_cast < std::string > ( code );
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";
    Error.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordErrorToErrorList( str );
                                                          
    // create context to write error in summary.map  
    str.erase(2,1);                                 
    S.str("");S.clear();                            
    S << ">>>> " +str +"    :    MAP error in file "; 
    S << __FILE__;                                   
    S << " on line ";                               
    S << __LINE__;                                  
    S << " occurred.";                              
    Msg.WriteErrorToOutputFile( S.str() );                     
  };//END try
};


/**
 * ============================================================================
 * setElementOptions
 *
 * @input : P            -- 
 * @input : T            -- 
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void MAP_OtherStateType_class::setElementOptions( Element_ptr       &P     , 
					    const std::string &T     , 
					    MAP_ErrStat       &Error , 
					    MAP_Message       &Msg ){
  EnumParser <ElementOptions> parser;
  std::string WORD = boost::to_upper_copy( T );

  switch( parser.parseElementOptions( WORD ) ){
    
  case PLOT :  // plot cable using matplotlib
    P->setPLOT_flag( true );
    break;
	
  case X_POS : // print node X position to output file
    P->setX_POS_flag( true );
    break;
	
  case Y_POS : // print node Y position to output file
    P->setY_POS_flag( true );
    break;

  case Z_POS : // print node Z position to output file
    P->setZ_POS_flag( true );
    break;

  case X_FORCE : // print node sum forces in the X direction to output file
    P->setX_FORCE_flag( true );
    break;

  case Y_FORCE : // print node sum forces in the Y direction to output file
    P->setY_FORCE_flag( true );
    break;

  case Z_FORCE : // print node sum forces in the Z direction to output file
    P->setZ_FORCE_flag( true );
    break;

  case LINE_TENSION : // print tension along the line to output file
    P->setLINE_TENSION_flag( true );
    break;
	
  case OMIT_CONTACT : // ignore cable/seabed contact
    P->setOMIT_CONTACT_flag( true );
    break;

  case LAY_LENGTH : // length of line touching the seafloor
    P->setLAY_LENGTH_flag( true );
    break;
  default :

    /** 
     * If the string is an empy space ' ', then avoid printing an error to the screen
     */
    if( boost::iequals( T , " " ) == false ){
      Error.set_error_key( MAP_WARNING );
      Msg.RecordWarningToWarningList("Could not interpret '" 
                       + T 
                       + "' as a valid element option in the MAP input file.");
    };// END if
    break;
  };//END switch
};


/**
 * ============================================================================
 * writeElementData
 *
 * @input : Msg          -- Error message status
 * ============================================================================
 */
void MAP_OtherStateType_class::writeElementData( MAP_Message &Msg ){
  // initlize the string we are writting to the output message
  std::string write = "";
  std::ostringstream S, first, last;
  S << std::fixed << std::setprecision(1);
  first << std::fixed << std::setprecision(1);
  last << std::fixed << std::setprecision(1);

  for ( unsigned int i=0 ; i<element.size() ; i++) {
    // Line 1: write the cable type to the Msg parameter
    Msg.WriteDataToOutputFile( "Element "
                        + boost::lexical_cast<std::string>( i+1 )
                        + " properties:\n"
                        +"-------------------------------------------------\n"
                        + "    Line type      : |  "
                        + element[i]->getName( )
                        + "\n");
	
    // Line 2: Write X position for the element
    // set the first and last node position 
    first << element[i]->getAnchorPosition( &Node::X );
    last << element[i]->getFairleadPosition( &Node::X );

    // this write the data to a message status string
    writeXYZData( "    X Position [m] : |  "       , 
                  first.str()                       , 
                  element[i]->getAnchorFlag( &Node::X ),
                  last.str()                        , 
                  element[i]->getFairleadFlag( &Node::X ),
                  Msg );

    // clear ostringstream for next use
    first.str(""); first.clear();
    last.str(""); last.clear();

    // Line 3: Write Y position for the element
    // set the first and last node position 
    first << element[i]->getAnchorPosition( &Node::Y );
    last << element[i]->getFairleadPosition( &Node::Y );

    // this write the data to a message status string
    writeXYZData( "    Y Position [m] : |  "       , 
                  first.str()                       , 
                  element[i]->getAnchorFlag( &Node::Y ),
                  last.str()                        , 
                  element[i]->getFairleadFlag( &Node::Y ),
                  Msg );

    // clear ostringstream for next use
    first.str(""); first.clear();
    last.str(""); last.clear();

    // Line 4: Write Z position for the element
    // set the first and last node position 
    first << element[i]->getAnchorPosition( &Node::Z );
    last << element[i]->getFairleadPosition( &Node::Z );

    // this write the data to a message status string
    writeXYZData( "    Z Position [m] : |  "       , 
                  first.str()                       , 
                  element[i]->getAnchorFlag( &Node::Z ),
                  last.str()                        , 
                  element[i]->getFairleadFlag( &Node::Z ),
                  Msg );

    // clear ostringstream for next use
    first.str(""); first.clear();
    last.str(""); last.clear();

    // Line 5: write the unstretched element length
    S << element[i]->getLu();
    if ( element[i]->getLuFlag() == false ) write +="(";// _TEXT_COLOR::BLUE;
    write += S.str();
	
    // Remove the blue text option
    if ( element[i]->getLuFlag() == false ) write +=")";//write += _TEXT_COLOR::END;

    Msg.WriteDataToOutputFile( "    Length     [m] : |  "
                        + write
                        + "\n");
    write.clear();
    S.str(""); S.clear();


    // Line 6: write the H horizontal element force
    S << element[i]->getH();
    if ( element[i]->getHFlag() == false ) write += "(";//_TEXT_COLOR::BLUE;
    write += S.str();
	
    // Remove the blue text option
    if ( element[i]->getHFlag() == false ) write += ")";//_TEXT_COLOR::END;

    Msg.WriteDataToOutputFile( "    H          [N] : |  "
                        + write
                        + "\n");
    write.clear();
    S.str(""); S.clear();	

    // Line 6: write V vertical element force
    S << element[i]->getV();
    if ( element[i]->getVFlag() == false ) write += "(";//_TEXT_COLOR::BLUE;
    write += S.str();
	
    // Remove the blue text option
    if ( element[i]->getVFlag() == false ) write += ")";//write += _TEXT_COLOR::END;

    Msg.WriteDataToOutputFile( "    V          [N] : |  "
                        + write
                        + "\n\n\n");
    write.clear();
    S.str(""); S.clear();	
  };//END for
};


/**
 * ============================================================================
 * writeXYZData
 *
 * @input : position     -- 
 * @input : tail         -- 
 * @input : tail_bool    -- 
 * @input : head         -- 
 * @input : head_bool    -- 
 * @input : Msg          -- Error message status
 * ============================================================================
 */
void MAP_OtherStateType_class::writeXYZData( const std::string &position , 
                                       const std::string &tail     ,
                                       const bool        tail_bool ,
                                       const std::string &head     ,
                                       const bool        head_bool ,
                                       MAP_Message       &Msg ){    
  std::string write = "";
  std::string temp = "";

  // Create the header for the value we are writting to the message string 
  write += position;

  // Add blue text if the tail variable is iterated
  if ( tail_bool == false ) temp += "(";

  // temp becomes tail. We don't want to modify tail directly 
  // because it is declared as a const std::string
  temp += tail;

  // add white spaces to the end of the temp string to occupy enough
  // space in the command prompt (to align columns together)
  // Remove the blue text option
  if ( tail_bool == false ) temp += ")";//_TEXT_COLOR::END;

  // this makes sure the do-while statement is valid
  assert( temp.size()<=10 );
  do {
    temp += " ";
  } while( temp.size()<11 );

  write += temp;     
  temp.clear();

  write += " -> ";

  // Add blue text if the head variable is iterated
  if ( head_bool == false ) temp += "(";//_TEXT_COLOR::BLUE;

  temp += head;

  // Remove the blue text option
  if ( head_bool == false ) temp += ")";//_TEXT_COLOR::END;

  // this makes sure the do-while statement is valid
  assert( temp.size()<=10 );

  do {
    temp += " ";
  } while( temp.size()<11 );

  write += temp;
  temp.clear();

  write += "\n";

  Msg.WriteDataToOutputFile( write );
};
//void MAP_OtherStateType_class::writeXYZData( const std::string &position  , 
//				       const std::string &tail      ,
//				       bool              tail_bool ,
//				       const std::string &head      ,
//				       bool              head_bool ,
//				       MAP_Message       &Msg ){    
//    std::string write = "";
//    std::string temp = "";
//
//    // Create the header for the value we are writting to the message string 
//    write += position;
//
//    // Add blue text if the tail variable is iterated
//    if ( tail_bool == false ) temp += "(";//_TEXT_COLOR::BLUE;
//    
//    /**
//     * temp becomes tail. We don't want to modify tail directly 
//     * because it is declared as a const std::string
//     */
//    temp += tail;
//    
//    /**
//     * add white spaces to the end of the temp string to occupy enough
//     * space in the command prompt (to align columns together)
//     * Remove the blue text option
//     */
//    if ( tail_bool == false ) temp += ")";//_TEXT_COLOR::END;
//
//    // this makes sure the do-while statement is valid
//    assert( temp.size()<=10 );
//    do{
//	temp += " ";
//    } while( temp.size()<11 );
//
//    write += temp;     
//    temp.clear();
//        
//    write += " -> ";
//
//    // Add blue text if the head variable is iterated
//    if ( head_bool == false ) temp += "(";//_TEXT_COLOR::BLUE;
//
//    temp += head;
//
//    // Remove the blue text option
//    if ( head_bool == false ) temp += ")";//_TEXT_COLOR::END;
//
//    // this makes sure the do-while statement is valid
//    assert( temp.size()<=10 );
//
//    do{
//	temp += " ";
//    } while( temp.size()<11 );
//    
//    write += temp;
//    temp.clear();
//    
//    write += "\n";
//
//    Msg.WriteDataToOutputFile( write );
//};


/**
 * ============================================================================
 * summary
 *
 * @output : 
 * ============================================================================
 */
std::string MAP_OtherStateType_class::summary( ) {

  // Now print the CableLibrary, Node and Element
  // variables to the Msg parameter

  MAP_Message T;

  T.WriteDataToOutputFile("\n");


  this->writeEnvironmentData          ( T );
  this->writeCableLibraryData         ( T );
  this->writeNodeData                 ( T );
  this->writeElementData              ( T );

  // Disable the ability to write the linearized stiffness matrix for now
  // until this is verified to work correctly
  //
  // If the numerics routine is initialized, run the 
  // stiffness matrix linearization
  if ( this->isNumericsUninitialized()==false ){ 
    this->writeLinearizedStiffnessMatrix( T );
  };// END if

  return T.GetCharacterString();
};


/**
 * ============================================================================
 * setCableLibraryReference
 *   
 * This function is always called to assign CableLibrary as input in the 
 * MAP_ParameterType_class object
 *
 * @input : T            -- 
 * @input : i            -- 
 * ============================================================================
 */
void MAP_OtherStateType_class::setCableLibraryReference( MAP_ParameterType_class &T , const int index ) { 
  // This pushes a reference of all CableLibrary variables into MAP_ParameterType_class
  T.pushVar( &property[index]->Diam );         // cable diameter
  T.pushVar( &property[index]->MassDenInAir ); // cable mass density
  T.pushVar( &property[index]->EA );           // cable stiffness
  T.pushVar( &property[index]->CB );           // cable/seabed friction coefficient
};


/**
 * ============================================================================
 * checkNodeReference
 *
 * @input : i            -- 
 * @input : ptr          -- 
 *
 * @output : 
 * ============================================================================
 */
int MAP_OtherStateType_class::checkNodeReference( const int index, VarType Node::* ptr ) const { 
  // @todo : remove assertion in next version of MAP
  assert( ((*node[index]).*ptr).reference_counter==1  ) ;
  return ((*node[index]).*ptr).reference_counter; 
};


/**
 * ============================================================================
 * getList
 *
 * This plots the "details" of the MAP_Other type in Python
 * ============================================================================
 */
std::string MAP_OtherStateType_class::getList(){
  std::cout << "attemping to get list\n";
  return this->OtherStateDataTypes.list();
}


/**
 * ============================================================================
 * checkElementReference
 *
 * @input : i            -- 
 * @input : prt          -- 
 *
 * @output : 
 * ============================================================================
 */
int MAP_OtherStateType_class::checkElementReference( const int index, VarType Element::* ptr ) const { 
  // @todo : remove assertion in next version of MAP
  assert( ((*element[index]).*ptr).reference_counter==1  ) ;
  return ((*element[index]).*ptr).reference_counter; 
};


/**
 * ============================================================================
 * addDepth
 *
 * @input : T            -- 
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void MAP_OtherStateType_class::addDepth( const std::string &T     ,
                                   MAP_ErrStat       &Error ,
                                   MAP_Message       &Msg) {
    
  if ( T == "" ) throw MAP_ERROR_7;

  // make sure the variable can be converted to a double
  try {
    this->depth = boost::lexical_cast<double>( T );
  } catch ( boost::bad_lexical_cast const& ) {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_10 );
    str += " ] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_10 );
        
    Error.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordErrorToErrorList( str );
                                                                                                     
    // create context to write error in summary.map                                          //          |
    str.erase(2,1);                                                                          //          |
    std::ostringstream S;                                                                    //          |
    S.str("");S.clear();                                                                     //          |
    S << ">>>> " +str +"    :    MAP error in file ";                                        //          |
    S << __FILE__;                                                                           //          |
    S << " on line ";                                                                        //          |
    S << __LINE__;                                                                           //          |
    S << " occurred.";                                                                       //          |
    Msg.WriteErrorToOutputFile( S.str() );                                                              //          |
  };//END try
};


/**
 * ============================================================================
 * addGravity 
 *
 * @input : T            -- 
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void MAP_OtherStateType_class::addGravity( const std::string &T     ,
                                     MAP_ErrStat       &Error ,
                                     MAP_Message       &Msg) {
    
  if ( T == "" ) throw MAP_ERROR_8;

  // make sure the variable can be converted to a double
  try {
    this->gravity = boost::lexical_cast<double>( T );
  }
  catch ( boost::bad_lexical_cast const& ) {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_11 );
    str += " ] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_11 );
        
    Error.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordErrorToErrorList( str );
                                                        
    // create context to write error in summary.map                                          //          |
    str.erase(2,1);                                                                          //          |
    std::ostringstream S;                                                                    //          |
    S.str("");S.clear();                                                                     //          |
    S << ">>>> " +str +"    :    MAP error in file ";                                        //          |
    S << __FILE__;                                                                           //          |
    S << " on line ";                                                                        //          |
    S << __LINE__;                                                                           //          |
    S << " occurred.";                                                                       //          |
    Msg.WriteErrorToOutputFile( S.str() );                                                              //          |
  };//END try
};


/**
 * ============================================================================
 * addSeaDensity
 *
 * @input : T            -- 
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void MAP_OtherStateType_class::addSeaDensity( const std::string &T     ,
                                        MAP_ErrStat       &Error ,
                                        MAP_Message       &Msg) {
    
  if ( T == "" ) throw MAP_ERROR_9;

  // make sure the variable can be converted to a double
  try {
    this->rho_sea = boost::lexical_cast<double>( T );
  }
  catch ( boost::bad_lexical_cast const& ) {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_12 );
    str += " ] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_12 );

    Error.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordErrorToErrorList( str );

    // create context to write error in summary.map                                          //          |
    str.erase(2,1);                                                                          //          |
    std::ostringstream S;                                                                    //          |
    S.str("");S.clear();                                                                     //          |
    S << ">>>> " +str +"    :    MAP error in file ";                                        //          |
    S << __FILE__;                                                                           //          |
    S << " on line ";                                                                        //          |
    S << __LINE__;                                                                           //          |
    S << " occurred.";                                                                       //          |
    Msg.WriteErrorToOutputFile( S.str() );                                                              //          |
  };//END try
};


/**
 * ============================================================================
 * setNodeReferenceToUserData
 *
 * Passes the address of a 'Connect' node to UserData so the Newton force 
 * balance equation can be access in the PETSc solving function
 *
 * @todo : a case may come up where one component of Newton's equations does not have
 *         to be solved.  When this comes us, we include a better method to count 
 *         the number of Newton equations solved.
 *
 * @input : i            -- 
 * ============================================================================
 */
void MAP_OtherStateType_class::setNodeReferenceToUserData( const int index ){ 
  Node* raw_Node = node[index].get();
  this->user_data.pushNodeBack( raw_Node ); 
};


/**
 * ============================================================================
 * setMessageReferenceToUserData
 *
 * @input :              -- 
 * ============================================================================
 */
void MAP_OtherStateType_class::setMessageReferenceToUserData( MAP_Message &Msg ) { 
  this->user_data.setMessage( Msg );
};


/**
 * ============================================================================
 * setErrorStatusReferenceToUserData
 *
 * @input :              -- 
 * ============================================================================
 */
void MAP_OtherStateType_class::setErrorStatusReferenceToUserData( MAP_ErrStat &error ) { 
  this->user_data.setErrorCode( error );
};


/**
 * ============================================================================
 * setElementReferenceToUserData
 *
 * Passes the adress of each element so that the continuous analytical catenary
 * equation can be accessed in the PETSc solve functions. 
 *
 * @todo : a case may come up where we may need only one (instead of two)
 *         catenary equations to solve the problem. This will come up when H, Lu
 *         is defined, but not V. In the future, write logic into the program so
 *         that the appropriate number of catenary equations are used to solve
 *         the problem. 
 *
 * @input : i            -- 
 * ============================================================================
 */
void MAP_OtherStateType_class::setElementReferenceToUserData( int i ) { 
  Element* raw_Element = element[i].get();
  //checkpoint();
  //std::cout << raw_Element << std::endl;
  this->user_data.pushElementBack( raw_Element );  
};


/**
 * ============================================================================
 * setMAP_ConstraintStateType_classReferenceToUserData
 *
 * Passes the address of the constraint states (the variables we are iterating/
 * solving)
 *
 * @input : T            -- 
 * ============================================================================
 */
void MAP_OtherStateType_class::setMAP_ConstraintStateType_classReferenceToUserData( MAP_ConstraintStateType_class &T ){ 
  this->user_data.initializeConstriant( T );
};


/**
 * ============================================================================
 * plot
 *
 * User option to plot the 3D cable profile with MatPlotLib
 *
 * @input : Error        -- Error code
 * @input : Msg          -- Error message status
 * ============================================================================
 */
void MAP_OtherStateType_class::plot( MAP_ErrStat &Error, MAP_Message &Msg ) {
    
  Msg.messageClean();

  std::vector <std::string> plot_string;

  std::string x     = "";//"pylab.plot(range(5))";
  std::string y     = "";
  std::string z     = "";
  std::string color = "";
  std::string out   = "";
  bool flag         = true;

  for( int i=0 ; i<this->getSizeOfElement() ; i++ ){
    // if the element has the option flag PLOT, then plot this element to the screen
    if ( this->getElementOptionFlag( i,&Element::PLOT_flag) ) {	    
      x     = "";
      y     = "";
      z     = "";
      color = "";
      out   = "";
      flag  = element[i]->elementPlotPoints(x,y,z,Msg);

      if ( flag != true ){
        Error.set_error_key( MAP_ERROR );
        Msg.RecordErrorToErrorList( "Plot aborted prematurely." );
      };// END if

      // get color
      if ( element[i]->getName( ) == "nylon" ){
        color = ",color=[.725,.725,.113],lw=3";  // dark yello
        //color = ",color='b',lw=3";
        //color = ",color=[.784,.568,.568],lw=3"; // pink
        //color = ",color=[.380,.709,.811],lw=3"; // blue
        //color = ",color=[0.015,0.749,1.0],lw=3"; // blue
        //color = ",color=[.898,.537,0],lw=3"; // orange
      }
      else if ( element[i]->getName( ) == "nylon2" ){
        //color = ",color=[.282,.487,.487],lw=3";
        color = ",color=[.917,.917,.674],lw=3"; // yello
        //color = ",color=[0,.7,1],lw=3";
      }
      else if (element[i]->getName( ) == "steel" ){
        color = ",color=[.725,.725,.113],lw=3";  // dark yello
        //color = ",color='b',lw=3";
        //color = ",color=[.725,.725,.113],lw=3";  // dark yello
        //color = ",color='r',lw=3";
      }
      else if (element[i]->getName( ) == "chainup" ){
        //color = ",color='k',lw=3";
      }
      else{ 
        color = ",'b'"; 
      }
    
      // add the string arrays into the Matplotlib plot command
      out += "ax.plot([";
      out += x;

      out += "],[";
      out += y;

      out += "],[";
      out += z;

      out += "]";
      out += color;//,color=[1,0,.35]
      out += ");";//,color=[1,0,.35]

//	    /**
//	     * =======  Plot black cirecly on line  ======       <--------------------------+
//	     */                                           		          //        |
//	                                                                          //        |
//	    std::vector <std::string> temp_string_x;                              //        |
//	    std::vector <std::string> temp_string_y;                              //        |
//	    std::vector <std::string> temp_string_z;                              //        |
//                                                                                  //        |
//	    boost::split(temp_string_x,x,boost::is_any_of(" , "));		  //        |
//	    boost::split(temp_string_y,y,boost::is_any_of(" , "));		  //        |
//	    boost::split(temp_string_z,z,boost::is_any_of(" , "));		  //        |
//										  //        |
//	    out += "ax.plot([";							  //        |
//	    out += temp_string_x[0] + " , " + temp_string_x[300] + "],[";	  //        |
//	    out += temp_string_y[0] + " , " + temp_string_y[300] + "],[";	  //        |
//	    out += temp_string_z[0] + " , " + temp_string_z[300] + "],";	  //        |
//	    out += "'ko')";							  //        |
//                                                                                  //        |
//	    //                                                   -----------------------------+ 
//	    //============== <END> ============================================================

      plot_string.push_back( out );
	    
    };//END if
  };

  if ( plot_string.size() != 0 ){    
    Py_Initialize();
    PyRun_SimpleString("import matplotlib as mpl");
//    PyRun_SimpleString("from matplotlib import rc");

//    PyRun_SimpleString("rc('text',usetex=True)");
//    PyRun_SimpleString("rc('font',family='serif')");

    //PyRun_SimpleString("rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})");
    //PyRun_SimpleString("rc('text', usetex=True)");
    PyRun_SimpleString("from mpl_toolkits.mplot3d import Axes3D");
    PyRun_SimpleString("import numpy as np");
    PyRun_SimpleString("import matplotlib.pyplot as plt");
//    PyRun_SimpleString("mpl.rcParams['legend.fontsize'] = 10");

    PyRun_SimpleString("fig = plt.figure()");
    PyRun_SimpleString("ax = Axes3D(fig)");

//    // Make background black
//    PyRun_SimpleString("ax.set_axis_bgcolor([0.13725,0.13725,0.13725])");

    // Now include the line points from Element::elementPlotPoints( ... )
    for ( unsigned int i=0 ; i<plot_string.size() ; i++ ){
      PyRun_SimpleString( plot_string[i].c_str() );
    };// END for

//    // Give axis limits
//    PyRun_SimpleString("ax.set_xlim3d(0, 400)");
//    PyRun_SimpleString("ax.set_ylim3d(-40, 40)");
//    PyRun_SimpleString("ax.set_zlim3d(-350, 100)");

//    // Make axis tick white
//    PyRun_SimpleString("for t in ax.w_xaxis.get_ticklines(): t.set_color('w')");
//    PyRun_SimpleString("for t in ax.w_yaxis.get_ticklines(): t.set_color('w')");
//    PyRun_SimpleString("for t in ax.w_zaxis.get_ticklines(): t.set_color('w')");
//    PyRun_SimpleString("for tick in ax.w_xaxis.get_major_ticks() + ax.w_yaxis.get_major_ticks() + ax.w_zaxis.get_major_ticks():\n for child in tick.get_children(): \n     child.set_color('w')");
//    PyRun_SimpleString("ax.set_xlabel(r'X Displacement [m]', color='w',fontsize=14)");
//    PyRun_SimpleString("ax.set_ylabel(r'Y Displacement [m]', color='w',fontsize=14)");
//    PyRun_SimpleString("ax.set_zlabel(r'Depth Below MSL [m]', color='w',fontsize=14)");


//    PyRun_SimpleString("ax.set_xlabel(r'X Displacement [m]', color='k',fontsize=14)");
//    PyRun_SimpleString("ax.set_ylabel(r'Y Displacement [m]', color='k',fontsize=14)");
//    PyRun_SimpleString("ax.set_zlabel(r'Depth Below MSL [m]', color='k',fontsize=14)");
    PyRun_SimpleString("ax.set_xlabel('X Displacement [m]', color='k',fontsize=14)");
    PyRun_SimpleString("ax.set_ylabel('Y Displacement [m]', color='k',fontsize=14)");
    PyRun_SimpleString("ax.set_zlabel('Depth Below MSL [m]', color='k',fontsize=14)");

//    // Plotting preferences and save the file
//    PyRun_SimpleString("ax.view_init(90, 0)");
//    PyRun_SimpleString("ax.view_init(30, 30)");
//    PyRun_SimpleString("fig.savefig('/media/sf_ISOPE_2013/MAP/Presentation/Figures/bridle2.pdf', transparent=True)");
	
    PyRun_SimpleString("plt.show()");
    //Py_Exit(0);
  };// END if
};


/**
 * ============================================================================
 * plot
 *
 * User option to plot the 3D cable profile with MatPlotLib
 *
 * @input : Error        -- Error code
 * @input : Msg          -- Error message status
 * ============================================================================
 */
std::vector <std::string> MAP_OtherStateType_class::plotString( MAP_ErrStat &Error, MAP_Message &Msg ) {

  MAP_Message T;

  std::vector <std::string> plot_string;

  std::string x     = "";
  std::string y     = "";
  std::string z     = "";
  std::string out   = "";
  bool flag         = true;

  for( int i=0 ; i<this->getSizeOfElement() ; i++ ){
    // if the element has the option flag PLOT, then plot this element to the screen
    if ( this->getElementOptionFlag( i,&Element::PLOT_flag) ) {	    
      x     = "";
      y     = "";
      z     = "";
      out   = "";
      flag  = element[i]->elementPlotPoints(x,y,z,Msg);

      if ( flag != true ){
        // @todo - add an error code that maps to the in the MAP_ERROR_n struct
        // and add a throw-catch statement. This looks like a hack job.
        Error.set_error_key( MAP_ERROR );
        Msg.RecordErrorToErrorList( "Plot aborted prematurely." );
      };// END if

      // add the string arrays into the Matplotlib plot command
      out += "[";
      out += x;
      plot_string.push_back( x );
//            std::vector <std::string> temp_string_x; 
//            boost::split(temp_string_x,x,boost::is_any_of(" , "));
//            for (int ii ; ii<temp_string_x.size() : ii++) std::cout << temp_string_x[ii] << std::endl;

      out += "],[";
      out += y;
      plot_string.push_back( y );

      out += "],[";
      out += z;
      out += "]";
      plot_string.push_back( z );

      //plot_string.push_back( out );
    };//END if
  };

  if ( plot_string.size() != 0 ){    
    // Now include the line points from Element::elementPlotPoints( ... )
    for ( unsigned int i=0 ; i<plot_string.size() ; i++ ){
      T.WriteDataToOutputFile( plot_string[i].c_str() );
    };// END for
  };// END if

  return plot_string;
};


/**
 * ============================================================================
 * setSolverOptions
 *
 * Reads the string inputs from the 'SOLVER OPTIONS' section of the MAP input
 * file and set the PETSc run-time options accordiningly.
 *
 * @input : T            -- 
 * @input : Msg          -- Error message status
 * ============================================================================
 */
void MAP_OtherStateType_class::setSolverOptions( const std::string &T, MAP_Message &Msg ){
  this->numeric_method->setNumericsOptionsString( T );
};


/**
 * ============================================================================
 * Solve
 *
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void MAP_OtherStateType_class::
Solve( MAP_ErrStat &Error , 
       MAP_Message &Msg )
{
  numeric_method->PetscSolve( *this, Error, Msg );
};


/**
 * ============================================================================
 * CheckConvergence
 *
 * Check the residuals to make sure the solution is actually reached. 
 *
 * @todo: should tol be a run-time option?
 * ============================================================================
 */
int MAP_OtherStateType_class::
CheckResidualConvergence( MAP_ErrStat        &Error , 
                          MAP_Message        &Msg )
{
  double tol = 1e-2;
  
  for ( unsigned int i=0 ; i<this->node.size() ; i++ ) {	
    // solve X direction Newton equation
    if (this->node[i]->getXNewtonEquationFlag()==true){      
      if (this->node[i]->f_x() >= tol) return 1;//throw MAP_ERROR_69;
      //std::cout << this->node[i]->f_x() << std::endl;
    };// END if
	
    // solve Y direction Newton equation
    if (this->node[i]->getYNewtonEquationFlag()==true){
      if (this->node[i]->f_y() >= tol) return 1;//throw MAP_ERROR_69;
      //std::cout << this->node[i]->f_y() << std::endl;
    };// END if

    // solve Z direction Newton equation
    if (this->node[i]->getZNewtonEquationFlag()==true){
      if (this->node[i]->f_z() >= tol) return 1;//throw MAP_ERROR_69;
      //std::cout << this->node[i]->f_z() << std::endl;
    };
  };//END for

  // For each element, the X,Y,Z catenary equation is solved.  
  for ( unsigned int i=0 ; i<this->element.size() ; i++ ) {
    if (this->element[i]->f_h() >= tol) return 1;//throw MAP_ERROR_69;
    //std::cout << this->element[i]->f_h() << std::endl;
    if (this->element[i]->f_v() >= tol) return 1;//throw MAP_ERROR_69;
    //std::cout << this->element[i]->f_v() << std::endl;
  };//END for
  //std::cout << "\n";

  return 0; // no error
};


/**
 * ============================================================================
 * initializeNumericSolver
 *
 * @input : Init         -- 
 * @input : Msg          -- Error message status
 * ============================================================================
 */
void MAP_OtherStateType_class::
initializeNumericSolver( MAP_InitInputType_class &Init  , 
                         MAP_ErrStat             &Error , 
                         MAP_Message             &Msg ) 
{
  numeric_method->InitializeSolver( *this , Init , Error, Msg );
};


/**
 * ============================================================================
 * cleanNumericSolver
 *
 * Terminates the PETSc solve routines and destroys all PETSc objects.  
 *
 * @input : Error        -- Error code
 * @input : Msg          -- Error message status
 * ============================================================================
 */
void MAP_OtherStateType_class::cleanNumericSolver( MAP_ErrStat &Error, MAP_Message &Msg ){
  numeric_method->End( Error, Msg );
};


/**
 * ============================================================================
 * getNumEquations
 *
 * @output : 
 * ============================================================================
 */
int MAP_OtherStateType_class::getNumEquations() {
  return this->num_equations;
};


/**
 * ============================================================================
 * initializeCableElement
 *
 * @input : index             -- 
 * ============================================================================
 */
void MAP_OtherStateType_class::initializeCableElement( const int   index , 
                                                 MAP_ErrStat &Error , 
                                                 MAP_Message &Msg ) { 
  this->element[index]->initializeElement( gravity, rho_sea , Error , Msg );
};


/**
 * ============================================================================
 * initializeCableNode >
 *
 * this sets the initial conditions for the Fix and Vessel nodes
 *
 * @input : index             -- node number 
 * ============================================================================
 */
void MAP_OtherStateType_class::initializeCableNode( const int index ){
  this->node[index]->initializeNode( );
};


/**
 * ============================================================================
 * initializeSumForce
 *
 * set Node.sum_FX, Node.sum_FY and Node.sum_FZ to zero
 * ============================================================================
 */
void MAP_OtherStateType_class::initializeSumForce( ) {
  for( unsigned int i=0 ; i<node.size() ; i++) {
    this->node[i]->setSumForceToZero();
  }
};


/**
 * ============================================================================
 * solve_sum_force_equation_in_direction_X >
 *
 * sets the 'solve_X_Newton_equation' flag in class Node
 *
 * @input : i            -- 
 * @input : T            -- 
 * ============================================================================
 */
void MAP_OtherStateType_class::setSolveSumForceEquationInDirectionX( const int index , bool T ){ 
  this->node[index]->setXNewtonEquationFlag( T ); 
}


/**
 * ============================================================================
 * solve_sum_force_equation_in_direction_Y
 *
 * sets the 'solve_Y_Newton_equation' flag in class Node
 *
 * @input : i            -- 
 * @input : T            -- 
 * ============================================================================
 */
void MAP_OtherStateType_class::setSolveSumForceEquationInDirectionY( const int index , bool T ){ 
  this->node[index]->setYNewtonEquationFlag( T ); 
}


/**
 * ============================================================================
 * solve_sum_force_equation_in_direction_Z
 *
 * sets the 'solve_Z_Newton_equation' flag in class Node
 *
 * @input : i            -- 
 * @input : T            -- 
 * ============================================================================
 */
void MAP_OtherStateType_class::setSolveSumForceEquationInDirectionZ( const int index , bool T ){ 
  this->node[index]->setZNewtonEquationFlag( T ); 
}


/**
 * ============================================================================
 * getSumForceEquationFlagInDirectionX
 *
 * @input : i            -- 
 *
 * @output : 
 * ============================================================================
 */
bool MAP_OtherStateType_class::getSumForceEquationFlagInDirectionX( const int index ){ 
  return this->node[index]->getXNewtonEquationFlag( ); 
}


/**
 * ============================================================================
 * getSumForceEquationFlagInDirectionY
 *
 * @input : i            -- 
 *
 * @output : 
 * ============================================================================
 */
bool MAP_OtherStateType_class::getSumForceEquationFlagInDirectionY( const int index ){ 
  return this->node[index]->getYNewtonEquationFlag( ); 
}


/**
 * ============================================================================
 * getSumForceEquationFlagInDirectionZ
 *
 * @input : i            -- 
 *
 * @output : 
 * ============================================================================
 */
bool MAP_OtherStateType_class::getSumForceEquationFlagInDirectionZ( const int index ){ 
  return this->node[index]->getZNewtonEquationFlag( ); 
}


/**
 * ============================================================================
 * getOutputStreamHeader
 *
 * @input : i            -- 
 * @input : Msg          -- Error message status
 * ============================================================================
 */
void MAP_OtherStateType_class::getOutputStreamHeader( const int index , MAP_Message &Msg ) {

  // print the simulation time
  if (index+1 == 1 ){
    this->output_string += "Time           ";
  };//END if

  if ( this->getElementOptionFlag( index , &Element::X_POS_flag) ){
    this->output_string += VarType::writeGenericVarType_name( index , element[index]->fairlead->X );
  };//END if

  if ( this->getElementOptionFlag( index , &Element::Y_POS_flag) ){
    this->output_string += VarType::writeGenericVarType_name( index , element[index]->fairlead->Y );
  };//END if

  if ( this->getElementOptionFlag( index , &Element::Z_POS_flag) ){
    this->output_string += VarType::writeGenericVarType_name( index , element[index]->fairlead->Z );
  };//END if

  if ( this->getElementOptionFlag( index , &Element::X_FORCE_flag) ){
    this->output_string += VarType::writeGenericVarType_name( index , element[index]->HX );
  };//END if

  if ( this->getElementOptionFlag( index , &Element::Y_FORCE_flag) ){
    this->output_string += VarType::writeGenericVarType_name( index , element[index]->HY );
  };//END if

  if ( this->getElementOptionFlag( index , &Element::Z_FORCE_flag) ){
    this->output_string += VarType::writeGenericVarType_name( index , element[index]->V );
  };//END if

  if ( this->getElementOptionFlag( index , &Element::LINE_TENSION_flag) ){
    this->output_string +=element[index]->getLineTensionHeader(index);
  };//END if

    // end the line
  if (index+1 >= this->getSizeOfElement() ){
    this->output_string += "\n";
  };//END if

};


/**
 * ============================================================================
 * getOutputStreamUnits
 *
 * @input : i            -- 
 * @input : Msg          -- Error message status
 * ============================================================================
 */
void MAP_OtherStateType_class::getOutputStreamUnits( const int index , MAP_Message &Msg ) {

  // print the simulation time
  if (index+1 == 1 ){
    this->output_string += "[s]            ";
  }

  if ( this->getElementOptionFlag( index , &Element::X_POS_flag) ){
    this->output_string += "[m]            " ;
  };//END if

  if ( this->getElementOptionFlag( index , &Element::Y_POS_flag) ){
    this->output_string += "[m]            " ;
  };//END if

  if ( this->getElementOptionFlag( index , &Element::Z_POS_flag) ){
    this->output_string += "[m]            " ;
  };//END if

  if ( this->getElementOptionFlag( index , &Element::X_FORCE_flag) ){
    this->output_string += "[N]            " ;
  };//END if

  if ( this->getElementOptionFlag( index , &Element::Y_FORCE_flag) ){
    this->output_string += "[N]            " ;
  };//END if

  if ( this->getElementOptionFlag( index , &Element::Z_FORCE_flag) ){
    this->output_string += "[N]            " ;
  };//END if

  if ( this->getElementOptionFlag( index , &Element::LINE_TENSION_flag) ){
    this->output_string += "[N]            " ;
    this->output_string += "[N]            " ;
    this->output_string += "[N]            " ;
    this->output_string += "[N]            " ;
    this->output_string += "[N]            " ;
    this->output_string += "[N]            " ;
    this->output_string += "[N]            " ;
    this->output_string += "[N]            " ;
    this->output_string += "[N]            " ;
    this->output_string += "[N]            " ;
  };//END if

    // end the line
  if (index+1 >= this->getSizeOfElement() ){
    this->output_string += "\n";
  }
};


/**
 * ============================================================================
 * getOutputStreamValue
 *  
 * @input : i            -- 
 * @input : Error        -- Error code 
 * ============================================================================
 */
void MAP_OtherStateType_class::getOutputStreamValue( const int index , const float time , MAP_Message &Msg ) {
// print the simulation time
//
// @todo : currently, only one time is printed. This will need to be changed
//         once MAP is included in a time-stepping routine. 

  if (index+1 == 1 ){
    std::ostringstream buff;
    buff << time;
    this->output_string += buff.str();
    do { 
      this->output_string += " ";
    } while( this->output_string.size()<15 ); 
  };// END if

  if ( this->getElementOptionFlag( index , &Element::X_POS_flag) ){
    this->output_string += VarType::writeGenericVarType_value( index , element[index]->fairlead->X );
  };//END if

  if ( this->getElementOptionFlag( index , &Element::Y_POS_flag) ){
    this->output_string += VarType::writeGenericVarType_value( index , element[index]->fairlead->Y );
  };//END if

  if ( this->getElementOptionFlag( index , &Element::Z_POS_flag) ){
    this->output_string += VarType::writeGenericVarType_value( index , element[index]->fairlead->Z );
  };//END if

  if ( this->getElementOptionFlag( index , &Element::X_FORCE_flag) ){
    this->output_string += VarType::writeGenericVarType_value( index , element[index]->HX );
  };//END if

  if ( this->getElementOptionFlag( index , &Element::Y_FORCE_flag) ){
    this->output_string += VarType::writeGenericVarType_value( index , element[index]->HY );
  };//END if

  if ( this->getElementOptionFlag( index , &Element::Z_FORCE_flag) ){
    this->output_string += VarType::writeGenericVarType_value( index , element[index]->V );
  };//END if

  if ( this->getElementOptionFlag( index , &Element::LINE_TENSION_flag) ){
    this->output_string += element[index]->getLineTension();
  };//END if

// end the line
  if (index+1 >= this->getSizeOfElement() ){
    this->output_string += "\n";
  };//END if
};




/**
 * ============================================================================
 * openOutputFile
 *
 * @todo : right now, MAP is hard coded to produce a 'map.out' file.  This is 
 *         probably not wanted; it does not allow users to name the output
 *         file whatever they want.  
 * ============================================================================
 */
void MAP_OtherStateType_class::openOutputFile( ){ 
  outputStream->open( "map.out" ); 
};


/**
 * ============================================================================
 * closeOutputFile
 *  
 *
 * ============================================================================
 */
void MAP_OtherStateType_class::closeOutputFile( ) { 
  *outputStream << "\n\nMAP successfully terminated";
  outputStream->close();
};


/**
 * ============================================================================
 * associateNodeVarTypeToNWTCType
 *
 * @input : i            -- 
 * @input : I            -- 
 * @input : P            -- 
 * @input : C            -- 
 * @input : O            -- 
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void MAP_OtherStateType_class::associateNodeVarTypeToNWTCType( const int               index ,
                                                         MAP_InputType_class           &I    ,
                                                         MAP_ParameterType_class       &P    , 
                                                         MAP_ConstraintStateType_class &C    ,
                                                         MAP_OutputType_class          &O    ,
                                                         MAP_ErrStat             &Err  ,
                                                         MAP_Message             &Msg){  

  /**
   * ==========   X Node VarType assignment  ===============     <-------------------------------------------------------------+
   * Assign data in the Node                                                                                       //          |
   */                                                                                                              //          |
  //          |
  try {                                                                                                            //          |
    // If X has a fixed value and it IS NOT attached to a vessel,                                                //          |
    // set it as a parameter                                                                                     //          |
    if ( this->nodeVarTypeBool( index , &Node::X ) == true && this->getNodeType( index ) != Vessel ) {           //          |
      // assign as a MAP_ParameterType_class                                                                         //          |
      this->setNodeReferenceToMAPType( index , P , &Node::X );                                                 //          |
    }                                                                                                            //          |
    // If X has a fixed value and it is attached to a vessel, set                                                //          |
    // it as an input                                                                                            //          |
    else if ( this->nodeVarTypeBool( index , &Node::X ) == true && this->getNodeType( index ) == Vessel ) {      //          |
      // assign as a MAP_InputType_class                                                                             //          |
      this->setNodeReferenceToMAPType( index , I , &Node::X );                                                 //          |
    }                                                                                                            //          |
    //          |
    // if X is variable (iterated), then it must be a constraint                                                 //          |
    else if ( this->nodeVarTypeBool( index , &Node::X ) == false ) {                                             //          |
      // ensure this node is not a 'Connect' or a 'Vessel'                                                     //          |
      if( this->getNodeType( index ) == Fix ){                                                                 //          |
        // @todo : remove assertion in the next release of MAP                                               //          |
        // once we verify things are working correctly                                                       //          |
        assert( this->getNodeType( index ) == Connect || this->getNodeType( index ) == Vessel );             //          |
        throw MAP_ERROR_27;                                                                                  //          |
      };//END if                                                                                               //          |
      //          |
      // assign as a MAP_ConstraintStateType_class                                                                   //          |
      this->setNodeReferenceToMAPType( index , C , &Node::X );                                                 //          |
      //          |
      // Increment the number of equations by                                                                  //          |
      // 1 only if the following conditions are                                                                //          |
      // met:                                                                                                  //          |
      //   -- X.is_fixed  = false                                                                              //          |
      //   -- FX.is_fixed = true                                                                               //          |
      //   -- Node.type   = Connect                                                                            //          |
      if ( this->nodeVarTypeBool( index, &Node::X )==false  && this->nodeVarTypeBool( index, &Node::FX )==true //          |
           && this->getNodeType( index ) == Connect ) {                                                        //          |
        // @todo : this is the world's worst way to increment                                                //          |
        // the number of equations we are solving. Let's fix                                                 //          |
        // this before we embarrass ourselves                                                                //          |
        this->incrementNumEquations();                                                                       //          |
        this->setSolveSumForceEquationInDirectionX( index , true );                                          //          |
      }                                                                                                        //          |
      else if ( this->getNodeType( index ) == Connect ) {                                                      //          |
        // The node displacement and applied for cannot                                                      //          |
        // be simultaneously defined. Here we issue a warning                                                //          |
        // message                                                                                           //          |
        throw MAP_ERROR_28;                                                                                  //          |
      };// END if                                                                                              //          |
    }                                                                                                            //          |
    else {                                                                                                       //          |
      // The node does not meet any initialization criteria                                                    //          |
      throw MAP_ERROR_29;                                                                                      //          |
    };//END if                                                                                                   //          |
  } catch ( MAP_ERROR_CODE &code ) {                                                                               //          |
    std::string str = "";                                                                                        //          |
    std::ostringstream S;                                                                                        //          |
    S << index+1;                                                                                                //          |
    str += boost::lexical_cast < std::string > ( code );                                                         //          |
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";                                         //          |
    //          |
    Err.set_error_key( MAP_ERROR );                                                                              //          |
    //          |
    // Here is the exception message                                                                             //          |
    Msg.RecordErrorToErrorList( str );                                                                                       //          |
    //          |
    // create context to write error in summary.map                                                              //          |
    str.erase(2,1);                                                                                              //          |
    S.str("");S.clear();                                                                                         //          |
    S << ">>>> " +str +"    :    MAP error in file ";                                                            //          |
    S << __FILE__;                                                                                               //          |
    S << " on line ";                                                                                            //          |
    S << __LINE__;                                                                                               //          |
    S << " occurred.";                                                                                           //          |
    Msg.WriteErrorToOutputFile( S.str() );                                                                                  //          |
  };//END try                                                                                                      //  --------+ 
    //============== <END> =======================================================================================================


    /**
     * ==========   Y Node VarType assignment  ===============     <-------------------------------------------------------------+
     * Assign data in the Node                                                                                       //          |
     */                                                                                                              //          |
                                                                                                                     //          |
  try {                                                                                                            //          |
    // If Y has a fixed value and it IS NOT attached to a vessel,                                                //          |
    // set it as a parameter                                                                                     //          |
    if ( this->nodeVarTypeBool( index , &Node::Y ) == true && this->getNodeType( index ) != Vessel ){            //          |
      // assign as a MAP_ParameterType_class                                                                         //          |
      this->setNodeReferenceToMAPType( index , P, &Node::Y );                                                  //          |
    }                                                                                                            //          |
    // If Y has a fixed value and it is attached to a vessel, set                                                //          |
    // it as an input                                                                                            //          |
    else if ( this->nodeVarTypeBool( index , &Node::Y ) == true && this->getNodeType( index ) == Vessel ){       //          |
      // assign as a MAP_InputType_class                                                                             //          |
      this->setNodeReferenceToMAPType( index , I , &Node::Y );                                                 //          |
    }                                                                                                            //          |
    //          |
    // if Y is variable (iterated), then it must be a constraint                                                 //          |
    else if ( this->nodeVarTypeBool( index , &Node::Y ) == false ) {                                             //          |
      // ensure this node is not a 'Connect' or a 'Vessel'                                                     //          |
      if( this->getNodeType( index ) == Fix ){                                                                 //          |
        // @todo : remove assertion in the next release of MAP                                               //          |
        // once we verify things are working correctly                                                       //          |
        assert( this->getNodeType( index ) == Connect || this->getNodeType( index ) == Vessel );             //          |
        throw MAP_ERROR_27;                                                                                  //          |
      };//END if                                                                                               //          |
      //          |
      // assign as a MAP_ConstraintStateType_class                                                                   //          |
      this->setNodeReferenceToMAPType( index , C , &Node::Y );                                                 //          |
      //          |
      // Increment the number of equations by                                                                  //          |
      // 1 only if the following conditions are                                                                //          |
      // met:                                                                                                  //          |
      //   -- Y.is_fixed  = false                                                                              //          |
      //   -- FY.is_fixed = true                                                                               //          |
      //   -- Node.type   = Connect                                                                            //          |
      if ( this->nodeVarTypeBool( index, &Node::Y )==false  && this->nodeVarTypeBool( index, &Node::FY )==true //          |
           && this->getNodeType( index ) == Connect ) {                                                        //          |
        // @todo : this is the world's worst way to increment                                                //          |
        // the number of equations we are solving. Let's fix                                                 //          |
        // this before we embarrass ourselves                                                                //          |
        this->incrementNumEquations();                                                                       //          |
        this->setSolveSumForceEquationInDirectionY( index , true );                                          //          |
      }                                                                                                        //          |
      else if ( this->getNodeType( index ) == Connect ) {                                                      //          |
        // The node displacement and applied for cannot                                                      //          |
        // be simultaneously defined. Here we issue a warning                                                //          |
        // message                                                                                           //          |
        throw MAP_ERROR_30;                                                                                  //          |
      };// END if                                                                                              //          |
    }                                                                                                            //          |
    else {                                                                                                       //          |
      // The node does not meet any initialization criteria                                                    //          |
      throw MAP_ERROR_31;                                                                                      //          |
    };//END if                                                                                                   //          |
  } catch ( MAP_ERROR_CODE &code ) {                                                                               //          |
    std::string str = "";                                                                                        //          |
    std::ostringstream S;                                                                                        //          |
    S << index+1;                                                                                                //          |
    str += boost::lexical_cast < std::string > ( code );                                                         //          |
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";                                         //          |
    //          |
    Err.set_error_key( MAP_ERROR );                                                                              //          |
    //          |
    // Here is the exception message                                                                             //          |
    Msg.RecordErrorToErrorList( str );                                                                                       //          |
    //          |
    // create context to write error in summary.map                                                              //          |
    str.erase(2,1);                                                                                              //          |
    S.str("");S.clear();                                                                                         //          |
    S << ">>>> " +str +"    :    MAP error in file ";                                                            //          |
    S << __FILE__;                                                                                               //          |
    S << " on line ";                                                                                            //          |
    S << __LINE__;                                                                                               //          |
    S << " occurred.";                                                                                           //          |
    Msg.WriteErrorToOutputFile( S.str() );                                                                                  //          |
  };//END try                                                                                                      //  --------+ 
    //============== <END> =======================================================================================================


    /**
     * ==========   Z Node VarType assignment  ===============     <-------------------------------------------------------------+
     * Assign data in the Node                                                                                       //          |
     */                                                                                                              //          |
                                                                                                                     //          |
  try {                                                                                                            //          |
    // If Z has a fixed value and it IS NOT attached to a vessel,                                                //          |
    // set it as a parameter                                                                                     //          |
    if ( this->nodeVarTypeBool( index , &Node::Z ) == true && this->getNodeType( index ) != Vessel ) {           //          |
      // assign as a MAP_ParameterType_class                                                                         //          |
      this->setNodeReferenceToMAPType( index , P, &Node::Z );                                                  //          |
    }                                                                                                            //          |
    // If Z has a fixed value and it is attached to a vessel, set                                                //          |
    // it as an input                                                                                            //          |
    else if ( this->nodeVarTypeBool( index , &Node::Z ) == true && this->getNodeType( index ) == Vessel ) {      //          |
      // assign as a MAP_InputType_class                                                                             //          |
      this->setNodeReferenceToMAPType( index , I , &Node::Z );                                                 //          |
    }                                                                                                            //          |
    //          |
    // if Z is variable (iterated), then it must be a constraint                                                 //          |
    else {                                                                                                       //          |
      // @todo : delete this in the next version of MAP after we verify this                                   //          |
      // logic is correct                                                                                      //          |
      assert( this->getNodeType( index ) == Connect );                                                         //          |
      //          |
      // ensure this node is a 'Connect'                                                                       //          |
      if( this->getNodeType( index ) != Connect ) {                                                            //          |
        // @todo : remove assertion in the next release of MAP                                               //          |
        // once we verify things are working correctly                                                       //          |
        //          |
        throw MAP_ERROR_27;                                                                                  //          |
      };//END if                                                                                               //          |
      //          |
      // assign as a MAP_ConstraintStateType_class                                                                   //          |
      this->setNodeReferenceToMAPType( index , C , &Node::Z );                                                 //          |
      //          |
      // Increment the number of equations by                                                                  //          |
      // 1 only if the following conditions are                                                                //          |
      // met:                                                                                                  //          |
      //   -- Z.is_fixed  = false                                                                              //          |
      //   -- FZ.is_fixed = true                                                                               //          |
      //   -- Node.type   = Connect                                                                            //          |
      if ( this->nodeVarTypeBool( index, &Node::Z )==false && this->nodeVarTypeBool( index, &Node::FZ )==true  //          |
           && this->getNodeType( index ) == Connect ) {                                                        //          |
        // @todo : this is the world's worst way to increment                                                //          |
        // the number of equations we are solving. Let's fix                                                 //          |
        // this before we embarrass ourselves                                                                //          |
        this->incrementNumEquations();                                                                       //          |
        this->setSolveSumForceEquationInDirectionZ( index , true );                                          //          |
      }                                                                                                        //          |
      else {                                                                                                   //          |
        // The node displacement and applied for cannot                                                      //          |
        // be simultaneously defined. Here we issue an error                                                 //          |
        // message                                                                                           //          |
        throw MAP_ERROR_32;                                                                                  //          |
      };// END if                                                                                              //          |
    };// END if                                                                                                  //          |
  } catch ( MAP_ERROR_CODE &code ) {                                                                               //          |
    std::string str = "";                                                                                        //          |
    std::ostringstream S;                                                                                        //          |
    S << index+1;                                                                                                //          |
    str += boost::lexical_cast < std::string > ( code );                                                         //          |
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";                                         //          |
    //          |
    Err.set_error_key( MAP_ERROR );                                                                              //          |
    //          |
    // Here is the exception message                                                                             //          |
    Msg.RecordErrorToErrorList( str );                                                                                       //          |
    //          |
    // create context to write error in summary.map                                                              //          |
    str.erase(2,1);                                                                                              //          |
    S.str("");S.clear();                                                                                         //          |
    S << ">>>> " +str +"    :    MAP error in file ";                                                            //          |
    S << __FILE__;                                                                                               //          |
    S << " on line ";                                                                                            //          |
    S << __LINE__;                                                                                               //          |
    S << " occurred.";                                                                                           //          |
    Msg.WriteErrorToOutputFile( S.str() );                                                                                  //          |
  };//END try                                                                                                      //  --------+ 
    //============== <END> =======================================================================================================


    /**
     * ==========   M Node VarType assignment  ===============     <-------------------------------------------------------------+
     * Assign data in the Node                                                                                       //          |
     */                                                                                                              //          |
                                                                                                                     //          |
  try {                                                                                                            //          |
    if ( this->nodeVarTypeBool( index , &Node::M ) == true ) {                                                   //          |
      // assign as a MAP_ParameterType_class                                                                         //          |
      this->setNodeReferenceToMAPType( index , P, &Node::M );                                                  //          |
    }                                                                                                            //          |
    else {                                                                                                       //          |
      // assign as a MAP_ConstraintStateType_class                                                                   //          |
      this->setNodeReferenceToMAPType( index , C , &Node::M );                                                 //          |
      //          |
      // Increment the number of equations by 1 only if the following conditions are met:                      //          |
      //   -- M.is_fixed  = false                                                                              //          |
      //   -- B.is_fixed  = true                                                                               //          |
      //   -- Z.is_fixe d = true                                                                               //          |
      //   -- FZ.is_fixed = true                                                                               //          |
      //   -- Node.type   = Connect                                                                            //          |
      if ( this->nodeVarTypeBool( index, &Node::M )==false && this->nodeVarTypeBool( index, &Node::Z )==true   //          |
           && this->nodeVarTypeBool( index, &Node::B )==true && this->nodeVarTypeBool( index, &Node::FZ )==true//          |
           && this->getNodeType( index ) == Connect ) {                                                        //          |
        this->incrementNumEquations();                                                                       //          |
      }                                                                                                        //          |
      else if ( this->getNodeType( index ) == Connect ) {                                                      //          |
        // The node displacement and applied for cannot be simultaneously defined. Here we                   //          |
        // issue a warning message                                                                           //          |
        throw MAP_ERROR_33;                                                                                  //          |
      };                                                                                                       //          |
      //          |
      if ( this->getNodeType( index ) != Fix ) {                                                               //          |
        throw MAP_ERROR_35;                                                                                  //          |
      };// END if                                                                                              //          |
    };// END if                                                                                                  //          |
  } catch ( MAP_ERROR_CODE &code ) {                                                                               //          |
    std::string str = "";                                                                                        //          |
    std::ostringstream S;                                                                                        //          |
    S << index+1;                                                                                                //          |
    str += boost::lexical_cast < std::string > ( code );                                                         //          |
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";                                         //          |
    //          |
    Err.set_error_key( MAP_ERROR );                                                                              //          |
    //          |
    // Here is the exception message                                                                             //          |
    Msg.RecordErrorToErrorList( str );                                                                                       //          |
    //          |
    //          |
    // create context to write error in summary.map                                                              //          |
    str.erase(2,1);                                                                                              //          |
    S.str("");S.clear();                                                                                         //          |
    S << ">>>> " +str +"    :    MAP error in file ";                                                            //          |
    S << __FILE__;                                                                                               //          |
    S << " on line ";                                                                                            //          |
    S << __LINE__;                                                                                               //          |
    S << " occurred.";                                                                                           //          |
    Msg.WriteErrorToOutputFile( S.str() );                                                                                  //          |
  };//END try                                                                                                      //  --------+ 
    //============== <END> =======================================================================================================


    /**
     * ==========   B Node VarType assignment  ===============     <-------------------------------------------------------------+
     * Assign data in the Node                                                                                       //          |
     */                                                                                                              //          |
                                                                                                                     //          |
  try {                                                                                                            //          |
    if ( this->nodeVarTypeBool( index , &Node::B ) == true ) {                                                   //          |
      // assign as a MAP_ParameterType_class                                                                         //          |
      this->setNodeReferenceToMAPType( index , P, &Node::B );                                                  //          |
    }                                                                                                            //          |
    else {                                                                                                       //          |
      // assign as a MAP_ConstraintStateType_class                                                                   //          |
      this->setNodeReferenceToMAPType( index , C , &Node::B );                                                 //          |
      //          |
      // Increment the number of equations by 1 only if the following conditions are met:                      //          |
      //   -- M.is_fixed  = true                                                                               //          |
      //   -- B.is_fixed  = false                                                                              //          |
      //   -- Z.is_fixed  = true                                                                               //          |
      //   -- FZ.is_fixed = true                                                                               //          |
      //   -- Node.type   = Connect                                                                            //          |
      if ( this->nodeVarTypeBool( index , &Node::B )==false && this->nodeVarTypeBool( index , &Node::Z )==true //          |
           && this->nodeVarTypeBool( index, &Node::M )==true && this->nodeVarTypeBool( index, &Node::FZ )==true//          |
           && this->getNodeType( index ) == Connect ) {                                                        //          |
        this->incrementNumEquations();                                                                       //          |
      }                                                                                                        //          |
      else if ( this->getNodeType( index ) == Connect ) {                                                      //          |
        // The node displacement and applied for cannot be simultaneously defined. Here we                   //          |
        // issue a warning message                                                                           //          |
        throw MAP_ERROR_34;                                                                                  //          |
      };                                                                                                       //          |
      //          |
      if ( this->getNodeType( index ) != Fix ) {                                                               //          |
        throw MAP_ERROR_35;                                                                                  //          |
      };// END if                                                                                              //          |
    };// END if                                                                                                  //          |
  } catch ( MAP_ERROR_CODE &code ) {                                                                               //          |
    std::string str = "";                                                                                        //          |
    std::ostringstream S;                                                                                        //          |
    S << index+1;                                                                                                //          |
    str += boost::lexical_cast < std::string > ( code );                                                         //          |
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";                                         //          |
    //          |
    Err.set_error_key( MAP_ERROR );                                                                              //          |
    //          |
    // Here is the exception message                                                                             //          |
    Msg.RecordErrorToErrorList( str );                                                                                       //          |
    //          |
    // create context to write error in summary.map                                                              //          |
    str.erase(2,1);                                                                                              //          |
    S.str("");S.clear();                                                                                         //          |
    S << ">>>> " +str +"    :    MAP error in file ";                                                            //          |
    S << __FILE__;                                                                                               //          |
    S << " on line ";                                                                                            //          |
    S << __LINE__;                                                                                               //          |
    S << " occurred.";                                                                                           //          |
    Msg.WriteErrorToOutputFile( S.str() );                                                                                  //          |
  };//END try                                                                                                      //  --------+ 
    //============== <END> =======================================================================================================


    /**
     * ===============  FX, FY and FZ rules  ======================================
     * 	   
     * This part is a little tricky. FX, FY and FZ can be set as iterated values in 
     * the MAP input file; however, they are not actually iterated (solved). If FX, 
     * FY or FZ is selected as an iterated value, it is a relfection of the sum forces 
     * applied to the node, H and V (once converted into the inertial/global frame).
     * 
     * The process below sets H and V at the element level as iterated variables.
     * 
     * @todo : incorporate M and B into this process
     *         once it is verified to work correctly
     * =============================================================================
     */


    /**
     * ==========   FX Node VarType assignment  ===============     <------------------------------------------------------------+
     * Assign data in the Node. FX can only be a parameter or constraint                                             //          |
     */                                                                                                              //          |
                                                                                                                     //          |
  try {                                                                                                            //          |
    if ( this->nodeVarTypeBool( index , &Node::FX ) == true ) {                                                  //          |
      // assign as a MAP_ParameterType_class                                                                         //

      this->setNodeReferenceToMAPType( index , P, &Node::FX );                                                 //          |
    }                                                                                                            //          |
    else {                                                                                                       //          |
      // @new
      //this->setNodeReferenceToMAPType( index , P , &Node::FX );                                                //          |
      this->setNodeReferenceToMAPType( index , this->OtherStateDataTypes , &Node::FX );

      // Increment the number of equations by 1 only if the following conditions are met:                      //          |
      //   -- FX.is_fixed = false                                                                              //          |
      //   -- X.is_fixed  = true                                                                               //          |
      //   -- Node.type   = Connect                                                                            //          |
      if ( this->nodeVarTypeBool( index, &Node::FX )==false && this->nodeVarTypeBool( index, &Node::X )==true  //          |
           && this->getNodeType( index ) == Connect ) {                                                        //          |
        this->incrementNumEquations();                                                                       //          |
      }                                                                                                        //          |
      else if ( this->getNodeType( index ) == Connect ) {                                                      //          |
        throw MAP_ERROR_36;                                                                                  //          |
      };                                                                                                       //          |
    };                                                                                                           //          |

    // if the node is attached to a vessel, then we must output                                                  //          |
    // the node force                                                                                            //          |
    if ( this->getNodeType( index ) == Vessel ){                                                                 //          |
      this->setNodeReferenceToMAPType( index , O , &Node::FX );                                                //          |
    };//END if                                                                                                   //          |
  } catch ( MAP_ERROR_CODE &code ) {                                                                               //          |
    std::string str = "";                                                                                        //          |
    std::ostringstream S;                                                                                        //          |
    S << index+1;                                                                                                //          |
    str += boost::lexical_cast < std::string > ( code );                                                         //          |
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";                                         //          |
    //          |
    Err.set_error_key( MAP_ERROR );                                                                              //          |
    //          |
    // Here is the exception message                                                                             //          |
    Msg.RecordErrorToErrorList( str );                                                                                       //          |
    //          |
    // create context to write error in summary.map                                                              //          |
    str.erase(2,1);                                                                                              //          |
    S.str("");S.clear();                                                                                         //          |
    S << ">>>> " +str +"    :    MAP error in file ";                                                            //          |
    S << __FILE__;                                                                                               //          |
    S << " on line ";                                                                                            //          |
    S << __LINE__;                                                                                               //          |
    S << " occurred.";                                                                                           //          |
    Msg.WriteErrorToOutputFile( S.str() );                                                                                  //          |
  };// END try                                                                                                     //----------+
    //============== <END> =======================================================================================================


    /**
     * ==========   FY Node VarType assignment  ===============     <------------------------------------------------------------+
     * Assign data in the Node. FY can only be a parameter or constraint                                             //          |
     */                                                                                                              //          |
                                                                                                                     //          |
  try {                                                                                                            //          |
    if ( this->nodeVarTypeBool( index , &Node::FY ) == true ) {                                                  //          |
      // assign as a MAP_ParameterType_class                                                                         //          |
      this->setNodeReferenceToMAPType( index , P, &Node::FY );                                                 //          |
    }                                                                                                            //          |
    else {                                                                                                       //          |
      // @new
      //this->setNodeReferenceToMAPType( index , P , &Node::FY );                                                //          |
      this->setNodeReferenceToMAPType( index , this->OtherStateDataTypes , &Node::FY );

      // Increment the number of equations by 1 only if the following conditions are met:                      //          |
      //   -- FY.is_fixed = false                                                                              //          |
      //   -- Y.is_fixed  = true                                                                               //          |
      //   -- Node.type   = Connect                                                                            //          |
      if ( this->nodeVarTypeBool( index, &Node::FY )==false && this->nodeVarTypeBool( index, &Node::Y )==true  //          |
           && this->getNodeType( index ) == Connect ) {                                                        //          |
        this->incrementNumEquations();                                                                       //          |
      }                                                                                                        //          |
      else if ( this->getNodeType( index ) == Connect ) {                                                      //          |
        throw MAP_ERROR_37;                                                                                  //          |
      };                                                                                                       //          |
    };                                                                                                           //          |
    //          |
    // if the node is attached to a vessel, then we must output                                                  //          |
    // the node force                                                                                            //          |
    if ( this->getNodeType( index ) == Vessel ){                                                                 //          |
      this->setNodeReferenceToMAPType( index , O , &Node::FY );                                                //          |
    };//END if                                                                                                   //          |
  } catch ( MAP_ERROR_CODE &code ) {                                                                               //          |
    std::string str = "";                                                                                        //          |
    std::ostringstream S;                                                                                        //          |
    S << index+1;                                                                                                //          |
    str += boost::lexical_cast < std::string > ( code );                                                         //          |
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";                                         //          |
    //          |
    Err.set_error_key( MAP_ERROR );                                                                              //          |
    //          |
    // Here is the exception message                                                                             //          |
    Msg.RecordErrorToErrorList( str );                                                                                       //          |
    //          |
    // create context to write error in summary.map                                                              //          |
    str.erase(2,1);                                                                                              //          |
    S.str("");S.clear();                                                                                         //          |
    S << ">>>> " +str +"    :    MAP error in file ";                                                            //          |
    S << __FILE__;                                                                                               //          |
    S << " on line ";                                                                                            //          |
    S << __LINE__;                                                                                               //          |
    S << " occurred.";                                                                                           //          |
    Msg.WriteErrorToOutputFile( S.str() );                                                                                  //          |
  };// END try                                                                                                     //----------+
    //============== <END> =======================================================================================================


    /**
     * ==========   FZ Node VarType assignment  ===============     <------------------------------------------------------------+
     * Assign data in the Node. FZ can only be a parameter or constraint                                             //          |
     */                                                                                                              //          |
                                                                                                                     //          |
  try {                                                                                                            //          |
    if ( this->nodeVarTypeBool( index , &Node::FZ ) == true ) {                                                  //          |
      // assign as a MAP_ParameterType_class                                                                         //          |
      this->setNodeReferenceToMAPType( index , P, &Node::FZ );                                                 //          |
    }                                                                                                            //          |
    else {                                                                                                       //          |
      // @new
      //this->setNodeReferenceToMAPType( index , P , &Node::FZ );                                                //          |
      this->setNodeReferenceToMAPType( index , this->OtherStateDataTypes , &Node::FZ );

      //          |
      // Increment the number of equations by 1 only if the following conditions are met:                      //          |
      //   -- FZ.is_fixed = false                                                                              //          |
      //   -- Z.is_fixed  = true                                                                               //          |
      //   -- B.is_fixed  = true                                                                               //          |
      //   -- M.is_fixed  = true                                                                               //          |
      //   -- Node.type   = Connect                                                                            //          |
      if ( this->nodeVarTypeBool( index, &Node::FZ )==false && this->nodeVarTypeBool( index, &Node::B )==true  //          |  
           && this->nodeVarTypeBool( index, &Node::M )==true && this->nodeVarTypeBool( index, &Node::Z )==true //          |  
           && this->getNodeType( index ) == Connect ) {                                                        //          |  
        this->incrementNumEquations();                                                                       //          |
      }                                                                                                        //          |
      else if ( this->getNodeType( index ) == Connect ) {                                                      //          |
        throw MAP_ERROR_38;                                                                                  //          |
      };                                                                                                       //          |
    };                                                                                                           //          |
    //          |
    // if the node is attached to a vessel, then we must output                                                  //          |
    // the node force                                                                                            //          |
    if ( this->getNodeType( index ) == Vessel ){                                                                 //          |
      this->setNodeReferenceToMAPType( index , O , &Node::FZ );                                                //          |
    };//END if                                                                                                   //          |
  } catch ( MAP_ERROR_CODE &code ) {                                                                               //          |
    std::string str = "";                                                                                        //          |
    std::ostringstream S;                                                                                        //          |
    S << index+1;                                                                                                //          |
    str += boost::lexical_cast < std::string > ( code );                                                         //          |
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";                                         //          |
    //          |
    Err.set_error_key( MAP_ERROR );                                                                              //          |
    //          |
    // Here is the exception message                                                                             //          |
    Msg.RecordErrorToErrorList( str );                                                                                       //          |
    //          |
    // create context to write error in summary.map                                                              //          |
    str.erase(2,1);                                                                                              //          |
    S.str("");S.clear();                                                                                         //          |
    S << ">>>> " +str +"    :    MAP error in file ";                                                            //          |
    S << __FILE__;                                                                                               //          |
    S << " on line ";                                                                                            //          |
    S << __LINE__;                                                                                               //          |
    S << " occurred.";                                                                                           //          |
    Msg.WriteErrorToOutputFile( S.str() );                                                                                  //          |
  };// END try                                                                                                     //----------+
    //============== <END> =======================================================================================================


    // final check to make sure we are not defining too many variables being solved 
    // (more equations than unknowns)
  int counter = 0;

  // Make sure Z direction doesn't have more unknown than equations
  try {
    if ( this->nodeVarTypeBool( index , &Node::M  ) == false ) counter++;
    if ( this->nodeVarTypeBool( index , &Node::Z  ) == false ) counter++;
    if ( this->nodeVarTypeBool( index , &Node::B  ) == false ) counter++;
    if ( this->nodeVarTypeBool( index , &Node::FZ ) == false ) counter++;

    if (counter > 1) {
      // The node displacement and applied for cannot        
      // be simultaneously defined. Here we issue a warning  
      // message                                             
      throw MAP_ERROR_39;
    };//END if
  } catch( MAP_ERROR_CODE &code ) {
    std::string str = "";
    std::ostringstream S;
    S << index+1;
    str += boost::lexical_cast < std::string > ( code );
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";

    Err.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordErrorToErrorList( str );
           
    // create context to write error in summary.map 
    str.erase(2,1);                                
    S.str("");S.clear();                            
    S << ">>>> " +str +"    :    MAP error in file ";
    S << __FILE__;                                  
    S << " on line ";                               
    S << __LINE__;                                  
    S << " occurred.";                              
    Msg.WriteErrorToOutputFile( S.str() );                     
  };// END try

    // Make sure X direction doesn't have more unknown than equations
  counter = 0;
  try {
    if ( this->nodeVarTypeBool( index , &Node::X  ) == false ) counter++;
    if ( this->nodeVarTypeBool( index , &Node::FX ) == false ) counter++;

    if (counter > 1) {
      // The node displacement and applied for cannot        
      // be simultaneously defined. Here we issue a warning  
      // message
      throw MAP_ERROR_40;
    };//END if
  } catch( MAP_ERROR_CODE &code ) {
    std::string str = "";
    std::ostringstream S;
    S << index+1;
    str += boost::lexical_cast < std::string > ( code );
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";

    Err.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordErrorToErrorList( str );

    // create context to write error in summary.map 
    str.erase(2,1);                                
    S.str("");S.clear();                            
    S << ">>>> " +str +"    :    MAP error in file ";
    S << __FILE__;                                  
    S << " on line ";                               
    S << __LINE__;                                  
    S << " occurred.";                              
    Msg.WriteErrorToOutputFile( S.str() );                     
  };// END try

    // Make sure Y direction doesn't have more unknown than equations
  counter = 0;
  try {
    if ( this->nodeVarTypeBool( index , &Node::Y  ) == false ) counter++;
    if ( this->nodeVarTypeBool( index , &Node::FY ) == false ) counter++;

    if (counter > 1){
      // The node displacement and applied for cannot be simultaneously defined. Here we issue a warning  
      // message            
      throw MAP_ERROR_41;
    };//END if 
  } catch( MAP_ERROR_CODE &code ) {
    std::string str = "";
    std::ostringstream S;
    S << index+1;
    str += boost::lexical_cast < std::string > ( code );
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";

    Err.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordErrorToErrorList( str );

    // create context to write error in summary.map 
    str.erase(2,1);                                
    S.str("");S.clear();                            
    S << ">>>> " +str +"    :    MAP error in file ";
    S << __FILE__;                                  
    S << " on line ";                               
    S << __LINE__;                                  
    S << " occurred.";                              
    Msg.WriteErrorToOutputFile( S.str() );                     
  };// END try
};


/**
 * ============================================================================
 * writeLinearizedStiffnessMatrix
 *
 * Computes the linearized stiffness matrix. This result will be writtent to
 * the MAP output file.
 *
 * The linearization model used centered finite differencing:
 *  
 *   -- [K]  =  [ F(X+epsilon) - F(X-epsilon) ]/(2*epsilon)
 *
 * @todo : In the future, the value for epsilon will be a run-time command. For
 *         now, it is set at compile time. 
 * @todo : In the cross-product term, check if the the element Hx, Hy and V 
 *         force should be used, or if the node force FX, FY and FZ should be 
 *         used. At the moment, FX, FY and HZ are being used. 
 * ============================================================================
 */
void MAP_OtherStateType_class::writeLinearizedStiffnessMatrix( MAP_Message &Msg ){

  std::string write        = "";
  std::string write_line_1 = "";
  std::string write_line_2 = "";
  std::string write_line_3 = "";
  std::string write_line_4 = "";
  std::string write_line_5 = "";
  std::string write_line_6 = "";


  bool        flg          = true;
  int         SIX          = 6; 
  double      epsilon      = 1e-3;
  double      delta        = 0.0;
  double      delta_2      = 0.0;
  double      force[6]     = { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 };
  double      Mx           = 0.0;
  double      My           = 0.0;
  double      Mz           = 0.0;

//    double X = 0.0;
//    double Y = 0.0;
//    double Z = 0.0;

  MAP_ErrStat        Err; 
  std::ostringstream S;


  S << std::fixed << std::setprecision(4);

  for ( unsigned int i=0 ; i<node.size() ; i++ ) {                              
    if( node[i]->type == Vessel ){			      
      /**								      
       * For Vessel fairleads: if FX or FY value does not begin with a '#', 
       * then we are not 						      
       */								      
      if ( node[i]->get_ptr( &Node::X )  == false ||
           node[i]->get_ptr( &Node::Y )  == false ||
           node[i]->get_ptr( &Node::Z )  == false ||
           node[i]->get_ptr( &Node::FX ) == true  ||
           node[i]->get_ptr( &Node::FY ) == true  ||
           node[i]->get_ptr( &Node::FZ ) == true ){
        flg = false;
      }	
    }								     
  } 

  /**
   * Perform the finite difference routine for each of the 6 DOFs
   *
   *   -- X
   *   -- Y
   *   -- Z
   *   -- phi   (roll)
   *   -- theta (pitch)
   *   -- psi   (yaw)
   */
  if (flg){
    for ( int jj=0 ; jj<SIX ; jj++ ) {

      /**
       * =======  Foward finite difference:  ======     <----------------------------------------------------------+
       *                                                                                               //          |
       *  F( X + epsilon )                                                                             //          |
       */                                                                                              //          |
      //          |
      for ( unsigned int i=0 ; i<node.size() ; i++ ) {                                                 //          |
        if ( node[i]->type == Vessel ) {								 //          |
          /**											 //          |
           * Store the original vessel displacement in some varaible. This will be restored	 //          |
           * once the new displacement+epsilon is solved. 					 //          |
           *											 //          |
           * perturb the vessel by epsilon							 //          |
           */											 //          |
          if ( jj == 0 ) { // X + epsilon								 //          |
            delta = node[i]->getPriorXValue( )+epsilon;					 //          |
            //delta = node[i]->getPtrValue( &Node::X )+epsilon;					 //          |
            node[i]->setPtrValue( &Node::X , delta );						 //          |
          }											 //          |
          else if ( jj == 1 ) { // Y + epsilon							 //          |
            delta = node[i]->getPriorYValue( )+epsilon;					 //          |
            //delta = node[i]->getPtrValue( &Node::Y )+epsilon;					 //          |
            node[i]->setPtrValue( &Node::Y , delta );						 //          |
          }											 //          |
          else if ( jj == 2 ) { // Z + epsilon							 //          |
            delta = node[i]->getPriorZValue( )+epsilon;					 //          |
            //delta = node[i]->getPtrValue( &Node::Z )+epsilon;					 //          |
            node[i]->setPtrValue( &Node::Z , delta );						 //          |
          }											 //          |
          else if ( jj == 3 ) { // phi + epsilon							 //          |
            /**											 //          |
             * -(phi*z) about j axis 								 //          |
             *  (phi*y) about k axis								 //          |
             */                             							 //          |
            delta   = node[i]->getPriorYValue( ) 						 //          |
              - node[i]->getPriorZValue( )*epsilon; 					 //          |
            delta_2 = node[i]->getPriorZValue( ) 						 //          |
              + node[i]->getPriorYValue( )*epsilon; 					 //          |
//			delta   = node[i]->getPtrValue( &Node::Y ) 						 //          |
//			    - node[i]->getPtrValue( &Node::Z )*epsilon; 					 //          |
//			delta_2 = node[i]->getPtrValue( &Node::Z ) 						 //          |
//			    + node[i]->getPtrValue( &Node::Y )*epsilon; 					 //          |
//			//          |
            node[i]->setPtrValue( &Node::Y , delta   ); 					 //          |
            node[i]->setPtrValue( &Node::Z , delta_2 ); 					 //          |
          }											 //          |
          else if ( jj == 4 ) { // theta + epsilon						 //          |
            /**											 //          |
             *  (theta*z) about i axis 								 //          |
             * -(theta*x) about k axis								 //          |
             */                             							 //          |
            delta   = node[i]->getPriorXValue( ) 						 //          |
              + node[i]->getPriorZValue( )*epsilon; 					 //          |
            delta_2 = node[i]->getPriorZValue( ) 						 //          |
              - node[i]->getPriorXValue( )*epsilon; 					 //          |
//			delta   = node[i]->getPtrValue( &Node::X ) 						 //          |
//			    + node[i]->getPtrValue( &Node::Z )*epsilon; 					 //          |
//			delta_2 = node[i]->getPtrValue( &Node::Z ) 						 //          |
//			    - node[i]->getPtrValue( &Node::X )*epsilon; 					 //          |
            //          |
            node[i]->setPtrValue( &Node::X , delta   ); 					 //          |
            node[i]->setPtrValue( &Node::Z , delta_2 ); 					 //          |
          }											 //          |
          else if ( jj == 5 ) { // psi + epsilon							 //          |
            /**											 //          |
             * -(psi*y) about i axis 								 //          |
             *  (psi*x) about j axis								 //          |
             */                             							 //          |
            delta   = node[i]->getPriorXValue( ) 						 //          |
              - node[i]->getPriorYValue( )*epsilon; 					 //          |
            delta_2 = node[i]->getPriorYValue( ) 						 //          |
              + node[i]->getPriorXValue( )*epsilon; 					 //          |                           							 //          |
//			delta   = node[i]->getPtrValue( &Node::X ) 						 //          |
//			    - node[i]->getPtrValue( &Node::Y )*epsilon; 					 //          |
//			delta_2 = node[i]->getPtrValue( &Node::Y ) 						 //          |
//			    + node[i]->getPtrValue( &Node::X )*epsilon; 					 //          |
            //          |
            node[i]->setPtrValue( &Node::X , delta   ); 					 //          |
            node[i]->setPtrValue( &Node::Y , delta_2 ); 					 //          |
          };//END if										 //          |
        };//END if                                                                                   //          |    
      };//END for                                            //----------------------------------------------------+ 
      //============== <END> =======================================================================================


      /**
       * restore original solution before epsilon displacement.
       *
       * 1) call numeric_method->PetscSolve( Error , Msg );
       * 
       * 2) reset the sum forces FX, FY and FZ in each node
       */
      Solve( Err, Msg );  
      for( int i=0 ; i<getSizeOfNode() ; i++) {
        initializeCableNode( i );
      };//END for


      /**
       * =======  Mooring force with Forward finite difference  ======     <---------------------------------------+
       *                                                                                               //          |
       * cross product to calculate moments :                                                          //          |
       *												 //          |
       *   | i   j   k  |      Y*Fz  - Z*Fy								 //          |
       *   | X   Y   Z  |  =   Z*Fx  - X*Fz								 //          |
       *   | Fx  Fy  Fz |      X*Fy  - Y*Fx								 //          |
       */												 //          |
													 //          |
      for ( unsigned int i=0 ; i<node.size() ; i++ ) {						 //          |
        if ( node[i]->type == Vessel ) {								 //          |
          force[0] += node[i]->getPtrValue( &Node::FX );  					 //          |
          force[1] += node[i]->getPtrValue( &Node::FY );						 //          |
          force[2] += node[i]->getPtrValue( &Node::FZ );						 //          |
          //          |
          //		Mx = ((node[i]->getPtrValue( &Node::FZ ))   *(node[i]->getPriorYValue( ))	   
          //		       - (node[i]->getPtrValue( &Node::FY ))*(node[i]->getPriorZValue( )));  
          //		My = ((node[i]->getPtrValue( &Node::FX ))   *(node[i]->getPriorZValue( ))	   
          //		       - (node[i]->getPtrValue( &Node::FZ ))*(node[i]->getPriorXValue( )));  
          //		Mz = ((node[i]->getPtrValue( &Node::FY ))   *(node[i]->getPriorXValue( ))	   
          //		       - (node[i]->getPtrValue( &Node::FX ))*(node[i]->getPriorYValue( )));  

          Mx = ((node[i]->getPtrValue( &Node::FZ ) ) *(node[i]->getPtrValue( &Node::Y ))	
                - (node[i]->getPtrValue( &Node::FY ))*(node[i]->getPtrValue( &Node::Z )));
          My = ((node[i]->getPtrValue( &Node::FX ))  *(node[i]->getPtrValue( &Node::Z ))	
                - (node[i]->getPtrValue( &Node::FZ ))*(node[i]->getPtrValue( &Node::X )));
          Mz = ((node[i]->getPtrValue( &Node::FY ))  *(node[i]->getPtrValue( &Node::X ))	
                - (node[i]->getPtrValue( &Node::FX ))*(node[i]->getPtrValue( &Node::Y )));	

          if ( jj == 0 || jj == 1 || jj == 2 ){
            Mx = ((node[i]->getPtrValue( &Node::FZ ))  *(node[i]->getPriorYValue( ))	   
                  - (node[i]->getPtrValue( &Node::FY ))*(node[i]->getPriorZValue( )));  
            My = ((node[i]->getPtrValue( &Node::FX ))  *(node[i]->getPriorZValue( ))	   
                  - (node[i]->getPtrValue( &Node::FZ ))*(node[i]->getPriorXValue( )));  
            Mz = ((node[i]->getPtrValue( &Node::FY ))  *(node[i]->getPriorXValue( ))	   
                  - (node[i]->getPtrValue( &Node::FX ))*(node[i]->getPriorYValue( )));  

            force[3] += Mx;
            force[4] += My;
            force[5] += Mz;
          }		
          else if ( jj == 3 ) {
            force[3] +=  Mx;
            force[4] +=  My*cos(epsilon) - Mz*sin(epsilon);	
            force[5] +=  My*sin(epsilon) + Mz*cos(epsilon);
          }
          else if ( jj == 4) {
            force[3] +=  Mx*cos(epsilon) + Mz*sin(epsilon);
            force[4] +=  My;	
            force[5] +=  -Mx*sin(epsilon) + Mz*cos(epsilon);
          } 
          else if ( jj == 5 ){
            force[3] +=  Mx*sin(epsilon) + My*cos(epsilon);	
            force[4] +=  Mx*cos(epsilon) - My*sin(epsilon);	
            force[5] +=  Mz;
          };//END if

          Mx = 0.0;
          My = 0.0;
          Mz = 0.0;
        };//END if											 //          |
      };//END for                                            //----------------------------------------------------+ 
      //============== <END> =======================================================================================
								
										
      /** 									
       * Restore the original position of the vessel nodes			
       * prior to perturbing it by epsilon					
       */									
      for ( unsigned int i=0 ; i<node.size() ; i++ ){				
        if ( node[i]->type == Vessel ){					
          node[i]->setPtrValue( &Node::X , node[i]->getPriorXValue() );	
          node[i]->setPtrValue( &Node::Y , node[i]->getPriorYValue() );	
          node[i]->setPtrValue( &Node::Z , node[i]->getPriorZValue() );	
        };//END if                                                           
      };//END for


      /**
       * =======  Backward finite difference:  ======     <--------------------------------------------------------+
       *                                                                                               //          |
       *  F( X - epsilon )                                                                             //          |
       */                                                                                              //          |
      //          |
      for ( unsigned int i=0 ; i<node.size() ; i++ ) {                                                 //          | 
        /**												 //          |
         * Only Vessel nodes are moded, because we want the stiffness to be relative		 //          |
         * to small vessel displacements								 //          |
         */												 //          |
        if ( node[i]->type == Vessel ) {								 //          |
          /**											 //          |
           * Store the original vessel displacement in some varaible. This will be restored	 //          |
           * once the new displacement+epsilon is solved. 					 //          |
           *											 //          |
           * perturb the vessel by epsilon							 //          |
           */											 //          |
          if ( jj == 0 ) { // X - epsilon								 //          |
            delta = node[i]->getPriorXValue( )-epsilon;					 //          |
            //delta = node[i]->getPtrValue( &Node::X )-epsilon;					 //          |
            node[i]->setPtrValue( &Node::X , delta );						 //          |
          }											 //          |
          else if ( jj == 1 ) { // Y - epsilon
            delta = node[i]->getPriorYValue( )-epsilon;					 //          |							 //          |
            //delta = node[i]->getPtrValue( &Node::Y )-epsilon;					 //          |
            node[i]->setPtrValue( &Node::Y , delta );						 //          |
          }											 //          |
          else if ( jj == 2 ) { // Z - epsilon							 //          |
            delta = node[i]->getPriorZValue( )-epsilon;					 //          |
            //delta = node[i]->getPtrValue( &Node::Z )-epsilon;					 //          |
            node[i]->setPtrValue( &Node::Z , delta );						 //          |
          }											 //          |
          else if ( jj == 3 ) { // phi  - epsilon							 //          |
            /**											 //          |
             *  (phi*z) about j axis 								 //          |
             * -(phi*y) about k axis								 //          |
             */                             							 //          |
            delta   = node[i]->getPriorYValue( ) 						 //          |
              + node[i]->getPriorZValue( )*epsilon; 					 //          |
            delta_2 = node[i]->getPriorZValue( ) 						 //          |
              - node[i]->getPriorYValue( )*epsilon; 					 //          |
//			delta   = node[i]->getPtrValue( &Node::Y ) 						 //          |
//			    + node[i]->getPtrValue( &Node::Z )*epsilon; 					 //          |
//			delta_2 = node[i]->getPtrValue( &Node::Z ) 						 //          |
//			    - node[i]->getPtrValue( &Node::Y )*epsilon; 					 //          |

            node[i]->setPtrValue( &Node::Y , delta   ); 					 //          |
            node[i]->setPtrValue( &Node::Z , delta_2 ); 					 //          |
          }											 //          |
          else if ( jj == 4 ) { // theta - epsilon						 //          |
            /**											 //          |
             * -(theta*z) about i axis 								 //          |
             *  (theta*x) about k axis								 //          |
             */                             							 //          |
            delta   = node[i]->getPriorXValue( ) 						 //          |
              - node[i]->getPriorZValue( )*epsilon; 					 //          |
            delta_2 = node[i]->getPriorZValue( ) 						 //          |
              + node[i]->getPriorXValue( )*epsilon; 					 //          |
//			delta   = node[i]->getPtrValue( &Node::X ) 						 //          |
//			    - node[i]->getPtrValue( &Node::Z )*epsilon; 					 //          |
//			delta_2 = node[i]->getPtrValue( &Node::Z ) 						 //          |
//			    + node[i]->getPtrValue( &Node::X )*epsilon; 					 //          |
            //          |
            node[i]->setPtrValue( &Node::X , delta   ); 					 //          |
            node[i]->setPtrValue( &Node::Z , delta_2 ); 					 //          |
          }											 //          |
          else if ( jj == 5 ) { // psi - epsilon							 //          |
            /**											 //          |
             *  (psi*y) about i axis 								 //          |
             * -(psi*x) about j axis								 //          |
             */  
            delta   = node[i]->getPriorXValue( ) 						 //          |
              + node[i]->getPriorYValue( )*epsilon; 					 //          |
            delta_2 = node[i]->getPriorYValue( ) 						 //          |
              - node[i]->getPriorXValue( )*epsilon; 					 //          |                           							 //          |
//			delta   = node[i]->getPtrValue( &Node::X ) 						 //          |
//			    + node[i]->getPtrValue( &Node::Y )*epsilon; 					 //          |
//			delta_2 = node[i]->getPtrValue( &Node::Y ) 						 //          |
//			    - node[i]->getPtrValue( &Node::X )*epsilon; 					 //          |
            //          |
            node[i]->setPtrValue( &Node::X , delta   ); 					 //          |
            node[i]->setPtrValue( &Node::Y , delta_2 ); 					 //          |
          };//END if										 //          |
        };//END if                                                                                   //          |     
      };//END for                                            //----------------------------------------------------+ 
      //============== <END> =======================================================================================
	
	
      /**
       * restore original solution before epsilon displacement.
       *
       * 1) call numeric_method->PetscSolve( Error , Msg );
       * 
       * 2) reset the sum forces FX, FY and FZ in each node
       */
      Solve( Err, Msg );  // call numeric_method->PetscSolve( Error , Msg );
      for( int i=0 ; i<getSizeOfNode() ; i++) {
        initializeCableNode( i );
      };//END for
	
	
      /**
       * =======  Mooring force with Backward finite difference  ======     <--------------------------------------+
       *                                                                                               //          |
       * cross product to calculate moments :                                                          //          |
       *												 //          |
       *   | i   j   k  |      Y*Fz  - Z*Fy								 //          |
       *   | X   Y   Z  |  =   Z*Fx  - X*Fz								 //          |
       *   | Fx  Fy  Fz |      X*Fy  - Y*Fx								 //          |
       */												 //          |
													 //          |
      for ( unsigned int i=0 ; i<node.size() ; i++ ) {						 //          |
        if ( node[i]->type == Vessel ) {								 //          |
          force[0] -= node[i]->getPtrValue( &Node::FX );						 //          |
          force[1] -= node[i]->getPtrValue( &Node::FY );						 //          |
          force[2] -= node[i]->getPtrValue( &Node::FZ );						 //          |
          //          |
          //		force[3] -= ((node[i]->getPtrValue( &Node::FZ ) ) *(node[i]->getPriorYValue( ))	         //          |
          //		     	     - (node[i]->getPtrValue( &Node::FY ))*(node[i]->getPriorZValue( )));	 //          |
          //		force[4] -= ((node[i]->getPtrValue( &Node::FX ))  *(node[i]->getPriorZValue( ))	         //          |
          //		     	     - (node[i]->getPtrValue( &Node::FZ ))*(node[i]->getPriorXValue( )));	 //          |
          //		force[5] -= ((node[i]->getPtrValue( &Node::FY ))  *(node[i]->getPriorXValue( ))	         //          |
          //			     - (node[i]->getPtrValue( &Node::FX ))*(node[i]->getPriorYValue( )));	 //          |

          Mx = ((node[i]->getPtrValue( &Node::FZ ) ) *(node[i]->getPtrValue( &Node::Y ))	
                - (node[i]->getPtrValue( &Node::FY ))*(node[i]->getPtrValue( &Node::Z )));
          My = ((node[i]->getPtrValue( &Node::FX ))  *(node[i]->getPtrValue( &Node::Z ))	
                - (node[i]->getPtrValue( &Node::FZ ))*(node[i]->getPtrValue( &Node::X )));
          Mz = ((node[i]->getPtrValue( &Node::FY ))  *(node[i]->getPtrValue( &Node::X ))	
                - (node[i]->getPtrValue( &Node::FX ))*(node[i]->getPtrValue( &Node::Y )));	

          if ( jj == 0 || jj == 1 || jj == 2 ){
            Mx = ((node[i]->getPtrValue( &Node::FZ ) ) *(node[i]->getPriorYValue( ))	
                  - (node[i]->getPtrValue( &Node::FY ))*(node[i]->getPriorZValue( )));	
            My = ((node[i]->getPtrValue( &Node::FX ))  *(node[i]->getPriorZValue( ))	
                  - (node[i]->getPtrValue( &Node::FZ ))*(node[i]->getPriorXValue( )));	
            Mz = ((node[i]->getPtrValue( &Node::FY ))  *(node[i]->getPriorXValue( ))
                  - (node[i]->getPtrValue( &Node::FX ))*(node[i]->getPriorYValue( )));	

            force[3] -= Mx;
            force[4] -= My;
            force[5] -= Mz;
          }
          else if ( jj == 3 ) {
            force[3] -=  Mx;
            force[4] -=  My*cos(epsilon) - Mz*sin(epsilon);	
            force[5] -=  My*sin(epsilon) + Mz*cos(epsilon);
          }
          else if ( jj == 4) {
            force[3] -=  Mx*cos(epsilon) + Mz*sin(epsilon);
            force[4] -=  My;	
            force[5] -= -Mx*sin(epsilon) + Mz*cos(epsilon);
          } 
          else if ( jj == 5 ){
            force[3] -=  Mx*sin(epsilon) + My*cos(epsilon);	
            force[4] -=  Mx*cos(epsilon) - My*sin(epsilon);	
            force[5] -=  Mz;
          };//END if

          Mx = 0.0;
          My = 0.0;
          Mz = 0.0;
        };// END if                                                                                  //          | 
      };//END for                                            //----------------------------------------------------+ 
      //============== <END> =======================================================================================

	
      /** 
       * Restore the original position of the vessel nodes
       * prior to perturbing it by epsilon
       */
      for ( unsigned int i=0 ; i<node.size() ; i++ ){
        if ( node[i]->type == Vessel ){
          node[i]->setPtrValue( &Node::X , node[i]->getPriorXValue() );
          node[i]->setPtrValue( &Node::Y , node[i]->getPriorYValue() );
          node[i]->setPtrValue( &Node::Z , node[i]->getPriorZValue() );
        };//END if
      };//END for


      /**
       * restore original solution before epsilon displacement.
       *
       * 1) call numeric_method->PetscSolve( Error , Msg );
       * 
       * 2) reset the sum forces FX, FY and FZ in each node
       */
      Solve( Err, Msg );  
      for( int i=0 ; i<getSizeOfNode() ; i++) {
        initializeCableNode( i );
      };//END for


      /**
       * ==========   X row of matrix [K]   ================     <------------------+
       *                                                            //             |               
       * Prints the X row linearized stiffness 		      //             |
       *							      //             |
       * To be sure the length of 'write_line_1' is not greater     //             |
       * than length we are allow to across, this assert statement  //             |
       * is raised.						      //             |
       *							      //             |
       * Fill spaces in the 'write' string with white spaces to     //             |
       * occupy area in the output string (to make the output	      //             |
       * print in alignment)					      //             |
       */							      //             |
								      //             |
      S << force[0]/(2*epsilon);				      //             |
      if ( force[0]/(2*epsilon) >= 0.0) write += " ";		      //             |
      write += S.str();					      //             |
      S.str("");S.clear(); 					      //             |
      //             |
      assert( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      do {							      //             |
        write += " ";					      //             |
      } while( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      write_line_1 += write;                                        //             |
      write.clear();                                                //   ----------+
      //============================================================================

	
      /**
       * ==========   Y row of matrix [K]   ================     <------------------+
       *                                                            //             |               
       * Prints the Y row linearized stiffness 		      //             |
       *							      //             |
       * To be sure the length of 'write_line_1' is not greater     //             |
       * than length we are allow to across, this assert statement  //             |
       * is raised.						      //             |
       *							      //             |
       * Fill spaces in the 'write' string with white spaces to     //             |
       * occupy area in the output string (to make the output	      //             |
       * print in alignment)					      //             |
       */							      //             |
								      //             |
      S << force[1]/(2*epsilon);                                    //             |
      if ( force[1]/(2*epsilon) >= 0.0) write += " ";		      //             |
      write += S.str();					      //             |
      S.str("");S.clear(); 					      //             |
      //             |
      assert( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      do {							      //             |
        write += " ";					      //             |
      } while( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      write_line_2 += write;                                        //             |
      write.clear();                                                //   ----------+
      //============================================================================

	
      /**
       * ==========   Z row of matrix [K]   ================     <------------------+
       *                                                            //             |               
       * Prints the Z row linearized stiffness 		      //             |
       *							      //             |
       * To be sure the length of 'write_line_1' is not greater     //             |
       * than length we are allow to across, this assert statement  //             |
       * is raised.						      //             |
       *							      //             |
       * Fill spaces in the 'write' string with white spaces to     //             |
       * occupy area in the output string (to make the output	      //             |
       * print in alignment)					      //             |
       */							      //             |
								      //             |
      S << force[2]/(2*epsilon);                                    //             |
      if ( force[2]/(2*epsilon) >= 0.0) write += " ";		      //             |
      write += S.str();					      //             |
      S.str("");S.clear(); 					      //             |
      //             |
      assert( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      do {							      //             |
        write += " ";					      //             |
      } while( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      write_line_3 += write;                                        //             |
      write.clear();                                                //   ----------+
      //============================================================================

	
      /**
       * ==========   Roll row of matrix [K]   ================     <--------------+
       *                                                            //             |               
       * Prints the Roll row linearized stiffness 		      //             |
       *							      //             |
       * To be sure the length of 'write_line_1' is not greater     //             |
       * than length we are allow to across, this assert statement  //             |
       * is raised.						      //             |
       *							      //             |
       * Fill spaces in the 'write' string with white spaces to     //             |
       * occupy area in the output string (to make the output	      //             |
       * print in alignment)					      //             |
       */							      //             |
								      //             |
      S << force[3]/(2*epsilon);                                    //             |
      if ( force[3]/(2*epsilon) >= 0.0) write += " ";		      //             |
      write += S.str();					      //             |
      S.str("");S.clear(); 					      //             |
      //             |
      assert( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      do {							      //             |
        write += " ";					      //             |
      } while( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      write_line_4 += write;                                        //             |
      write.clear();                                                //   ----------+
      //============================================================================

	
      /**
       * ==========   Pitch row of matrix [K]   ================     <-------------+
       *                                                            //             |               
       * Prints the Pitch row linearized stiffness 		      //             |
       *							      //             |
       * To be sure the length of 'write_line_1' is not greater     //             |
       * than length we are allow to across, this assert statement  //             |
       * is raised.						      //             |
       *							      //             |
       * Fill spaces in the 'write' string with white spaces to     //             |
       * occupy area in the output string (to make the output	      //             |
       * print in alignment)					      //             |
       */							      //             |
								      //             |
      S << force[4]/(2*epsilon);                                    //             |
      if ( force[4]/(2*epsilon) >= 0.0) write += " ";		      //             |
      write += S.str();					      //             |
      S.str("");S.clear(); 					      //             |
      //             |
      assert( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      do {							      //             |
        write += " ";					      //             |
      } while( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      write_line_5 += write;                                        //             |
      write.clear();                                                //   ----------+
      //============================================================================

	
      /**
       * ==========   Yaw row of matrix [K]   ================     <---------------+
       *                                                            //             |               
       * Prints the Yaw row linearized stiffness 		      //             |
       *							      //             |
       * To be sure the length of 'write_line_1' is not greater     //             |
       * than length we are allow to across, this assert statement  //             |
       * is raised.						      //             |
       *							      //             |
       * Fill spaces in the 'write' string with white spaces to     //             |
       * occupy area in the output string (to make the output	      //             |
       * print in alignment)					      //             |
       */							      //             |
								      //             |
      S << force[5]/(2*epsilon);                                    //             |
      if ( force[5]/(2*epsilon) >= 0.0) write += " ";		      //             |
      write += S.str();					      //             |
      S.str("");S.clear(); 					      //             |
      //             |
      assert( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      do {							      //             |
        write += " ";					      //             |
      } while( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      write_line_6 += write;                                        //             |
      write.clear();                                                //   ----------+
      //============================================================================


      for( int i=0 ; i<SIX ; i++) {
        force[i] = 0.0;
      };//END for
    };//END for

    Msg.WriteDataToOutputFile( "Linearized Stiffness Matrix\n" );
    Msg.WriteDataToOutputFile( "--------------------------------------------------------------------------" );
    Msg.WriteDataToOutputFile( "--------------------------------------------------------------------------\n" );
    Msg.WriteDataToOutputFile( "    " + write_line_1 );
    Msg.WriteDataToOutputFile( "\n" );
    Msg.WriteDataToOutputFile( "    " + write_line_2 );
    Msg.WriteDataToOutputFile( "\n" );
    Msg.WriteDataToOutputFile( "    " + write_line_3 );
    Msg.WriteDataToOutputFile( "\n" );
    Msg.WriteDataToOutputFile( "    " + write_line_4 );
    Msg.WriteDataToOutputFile( "\n" );
    Msg.WriteDataToOutputFile( "    " + write_line_5 );
    Msg.WriteDataToOutputFile( "\n" );
    Msg.WriteDataToOutputFile( "    " + write_line_6 );
    Msg.WriteDataToOutputFile( "\n\n" );
  }
  else {
    Msg.WriteDataToOutputFile( "Linearized Stiffness Matrix\n" );
    Msg.WriteDataToOutputFile( "--------------------------------------------------------------------------" );
    Msg.WriteDataToOutputFile( "--------------------------------------------------------------------------\n" );
    Msg.WriteDataToOutputFile( "    Stiffness matrix could not be computed with given\n" );
    Msg.WriteDataToOutputFile( "    input file options\n" );
  }
};
