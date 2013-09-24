 /*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   MAP_OtherStateType.cpp
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


#include "MAP_OtherStateType.h"


// ============================================================================
// addCableLibrary
//
// Adds a cable type (such as steel, nylon, vectran) to the cable library
//
// @input : T[0]         -- name of variable
// @input : T[1]         -- diameter
// @input : T[2]         -- density of material
// @input : T[3]         -- Young's modulus
// @input : T[4]         -- friction coef. between seabed and cable 
// @input : Msg          -- Error message status
// @input : Error        -- Error code
// ============================================================================
void MAP_OtherStateType_class::
addCableLibrary( const std::vector<std::string> &T     ,
                 const int                      index  ,
                 MAP_ErrStat_class              &Error ,
                 MAP_Message_class              &Msg   ) 
{    
  CableLibrary_ptr prop_ptr( new CableLibrary );  // create new memory on the heap (dynamically allocated)
  prop_ptr->label = T[0];                         // First give the cable a name, such as steel or nylon

  /*
      Set cable diameter properties       <------------------------------------------------------+
  */
  prop_ptr->Diam.name = "Diam"; 

  if ( boost::starts_with( T[1] , "#" ) ) { // if the first elemt in diameter starts with '#'
    std::string str = "";
    MAPSetUniversalErrorStat( MAP_WARNING_2, str, Error, Msg );

    /*
        Now cycle throw the string and write it to the This parameter, but ignore the 
        first character (which will be a '#' due to if-statement condition above.
    */
    std::string This = "";                                                 
    for ( unsigned int i=1 ; i<T[1].size() ; i++) {                        
      This.push_back( T[1][i] );                                           
    }                                                                      

    try { // make sure the variable can be converted to a double           
      prop_ptr->Diam.value = boost::lexical_cast<double>( This );          
      prop_ptr->Diam.index = index;
    } catch ( boost::bad_lexical_cast& ) {                                 
      std::string str = " ";                                               
      str += boost::lexical_cast < std::string > ( MAP_ERROR_6 );          
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_6 );          

      Error.set_error_key( MAP_ERROR );                                    

      // Here is the exception message                                     
      Msg.RecordToErrorList( str );                                        

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
    }
  } else {
    try { // make sure the variable can be converted to a double           
      prop_ptr->Diam.value = boost::lexical_cast<double>( T[1] );          
      prop_ptr->Diam.index = index;
    } catch ( const boost::bad_lexical_cast& ) {                           
      std::string str = " ";                                               
      str += boost::lexical_cast < std::string > ( MAP_ERROR_6 );          
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_6 );          

      Error.set_error_key( MAP_ERROR );                                    

      // Here is the exception message                                     
      Msg.RecordToErrorList( str );                                        

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
    }
  }

  /*
      Set cable density properties       <-------------------------------------------------------+
  */
  prop_ptr->MassDenInAir.name = "MassDenInAir";  // set density properties  

  if ( boost::starts_with( T[2] , "#" ) ) {  // if the first character in MassDenInAir starts with '#' 
    std::string str = "";
    MAPSetUniversalErrorStat( MAP_WARNING_3, str, Error, Msg );

    // Now cycle throw the string and write it to the This parameter,   
    // but ignore the first character (which will be a '#' due to       
    // if-statement condition above                                     
    std::string This = "";                                              
    for ( unsigned int i=1 ; i<T[2].size() ; i++) {                     
      This.push_back( T[2][i] );                                        
    }

    try { // make sure the variable can be converted to a double        
      prop_ptr->MassDenInAir.value = boost::lexical_cast<double>( This );
      prop_ptr->MassDenInAir.index = index;
    } catch ( boost::bad_lexical_cast& ) {                              
      std::string str = "";                                             
      str += boost::lexical_cast < std::string > ( MAP_ERROR_17 );      
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_17 );      

      Error.set_error_key( MAP_ERROR );                                 

      // Here is the exception message                                  
      Msg.RecordToErrorList( str );                                     

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
    }
  } else {                                                                            
    try { // make sure the variable can be converted to a double                    
      prop_ptr->MassDenInAir.value = boost::lexical_cast<double>( T[2] );           
      prop_ptr->MassDenInAir.index = index;                                         
    } catch ( boost::bad_lexical_cast& ) {                                          
      std::string str = "";                                                         
      str += boost::lexical_cast < std::string > ( MAP_ERROR_17 );                  
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_17 );                  

      Error.set_error_key( MAP_ERROR );                                             

      // Here is the exception message                                              
      Msg.RecordToErrorList( str );                                                 

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
    }
  }

  /*
      Set cable Young's modulus properties       <-----------------------------------------------+
  */
  prop_ptr->EA.name = "EA";     

  if ( boost::starts_with( T[3] , "#" ) ) { // if the first character in 'E' starts with '#'   
    std::string str = "";
    MAPSetUniversalErrorStat( MAP_WARNING_4, str, Error, Msg );

    // Now cycle throw the string and write it to the This parameter,       
    // but ignore the first character (which will be a '#' due to           
    // if-statement condition above                                         
    std::string This = "";                                                  
    for ( unsigned int i=1 ; i<T[3].size() ; i++) {                         
      This.push_back( T[3][i] );                                            
    };// END for                                                            

    try { // make sure the variable can be converted to a double            
      prop_ptr->EA.value = boost::lexical_cast<double>( This );             
      prop_ptr->EA.index = index;
    } catch ( boost::bad_lexical_cast& ) {                                  
      std::string str = "";                                                 
      str += boost::lexical_cast < std::string > ( MAP_ERROR_18 );          
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_18 );          

      Error.set_error_key( MAP_ERROR );                                     

      // Here is the exception message                                      
      Msg.RecordToErrorList( str );                                         

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
    }
  }                                                                         
  else {                                                                    
    try { // make sure the variable can be converted to a double            
      prop_ptr->EA.value = boost::lexical_cast<double>( T[3] );             
      prop_ptr->EA.index = index;
    } catch ( boost::bad_lexical_cast& ) {                                  
      std::string str = "";                                                 
      str += boost::lexical_cast < std::string > ( MAP_ERROR_18 );          
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_18 );          

      Error.set_error_key( MAP_ERROR );                                     

      // Here is the exception message                                      
      Msg.RecordToErrorList( str );                                         

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
    }
  }

  /*
      Set cable/sea bed friction coefficients       <--------------------------------------------+
  */
  prop_ptr->CB.name = "CB";                                                                       

  if ( boost::starts_with( T[4] , "#" ) ) { // if the first character in 'CB starts with '#'      
    std::string str = "";
    MAPSetUniversalErrorStat( MAP_WARNING_5, str, Error, Msg );

    // Now cycle throw the string and write it to the This parameter,       
    // but ignore the first character (which will be a '#' due to           
    // if-statement condition above                                         
    std::string This = "";                                                  
    for ( unsigned int i=1 ; i<T[4].size() ; i++) {                         
      This.push_back( T[4][i] );                                            
    }

    try { // make sure the variable can be converted to a double            
      prop_ptr->CB.value = boost::lexical_cast<double>( This );             
      prop_ptr->CB.index = index;
    } catch ( boost::bad_lexical_cast& ) {                                  
      std::string str = "";                                                 
      str += boost::lexical_cast < std::string > ( MAP_ERROR_19 );          
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_19 );          

      Error.set_error_key( MAP_ERROR );                                     

      // Here is the exception message                                      
      Msg.RecordToErrorList( str );                                         

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
    }
  } else {                                                                    
    try { // make sure the variable can be converted to a double            
      prop_ptr->CB.value = boost::lexical_cast<double>( T[4] );             
      prop_ptr->CB.index = index;
    } catch ( boost::bad_lexical_cast& ) {                                  
      std::string str = "";                                                 
      str += boost::lexical_cast < std::string > ( MAP_ERROR_19 );          
      str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_19 );          

      Error.set_error_key( MAP_ERROR );                                     

      // Here is the exception message                                      
      Msg.RecordToErrorList( str );                                         

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
    }
  }

  property.push_back( prop_ptr ); // push the new object we created into our vector container CableLibrary
};


// ============================================================================
//
// ============================================================================
void MAP_OtherStateType_class::
SetFastCouplingFlag( const bool flag ) 
{    
  is_coupled_to_FAST = flag;
}


// ============================================================================
//
// ============================================================================
bool MAP_OtherStateType_class::
GetFastCouplingFlag( ) const 
{    
  return is_coupled_to_FAST;
}


/**
 * MAPCALL_GetWriteOutput
 *
 * passes array of floats so Fortran code can write out the output independantly of FAST
 */
void MAP_OtherStateType_class::
MAPCALL_GetWriteOutput( float * arr , const int len ) {
  int count = 0;
  
  for ( unsigned int i=0 ; i<element.size() ; i++ ) {
    if ( this->getElementOptionFlag( i , &Element::X_POS_flag) ) {
      arr[count] = static_cast<float>(element[i]->fairlead->X.value);
      count++;
    }
    
    if ( this->getElementOptionFlag( i , &Element::Y_POS_flag) ) {
      arr[count] = static_cast<float>(element[i]->fairlead->Y.value);
      count++;
    }
    
    if ( this->getElementOptionFlag( i , &Element::Z_POS_flag) ) {
      arr[count] = static_cast<float>(element[i]->fairlead->Z.value);
      count++;
    }
    
    if ( this->getElementOptionFlag( i , &Element::X_FORCE_flag) ) {
      arr[count] = static_cast<float>(element[i]->fairlead->FX.value);
      count++;
    }

    if ( this->getElementOptionFlag( i , &Element::Y_FORCE_flag) ) {
      arr[count] = static_cast<float>(element[i]->fairlead->FY.value);
      count++;
    }
    
    if ( this->getElementOptionFlag( i , &Element::Z_FORCE_flag) ) {
      arr[count] = static_cast<float>(element[i]->fairlead->FZ.value);
      count++;
    }
    
    if ( this->getElementOptionFlag( i , &Element::LAY_LENGTH_flag) ) {
      arr[count] = element[i]->GetLayLengthValue();
      count++;
    }
    
    if ( this->getElementOptionFlag( i , &Element::LINE_TENSION_flag) ) {
      float tempArr[10];
      element[i]->GetLineTensionValues( tempArr );

      arr[count  ] = tempArr[0];
      arr[count+1] = tempArr[1];
      arr[count+2] = tempArr[2];
      arr[count+3] = tempArr[3];
      arr[count+4] = tempArr[4];
      arr[count+5] = tempArr[5];
      arr[count+6] = tempArr[6];
      arr[count+7] = tempArr[7];
      arr[count+8] = tempArr[8];
      arr[count+9] = tempArr[9];
      count+=10;
    }

    if ( this->getElementOptionFlag( i , &Element::FAIR_TENSION_flag) ) {
      arr[count] = element[i]->GetFairleadTensionMagnitude();
      count++;
    }

    if ( this->getElementOptionFlag( i , &Element::ANCH_TENSION_flag) ) {
      arr[count] = element[i]->GetAnchorTensionMagnitude();
      count++;
    }
  }
  
  assert( count == len );
};


// ============================================================================
// writeCableLibraryData
//
// write all data from the CableLibrary class to the MAP_Message_class parameter
//
// @input : Msg          -- Error message status
// ============================================================================
void MAP_OtherStateType_class::
writeCableLibraryData( MAP_Message_class &Msg) 
{    
  for ( unsigned int i=0 ; i<property.size() ; i++){

    std::ostringstream S;
    S << std::fixed << std::setprecision(3);
	
    // Line 1 : 
    // write the cable type to the Msg parameter
    Msg.WriteDataToOutputFile("    Cable Type          : ");
    Msg.WriteDataToOutputFile( property[i]->label );
    Msg.WriteDataToOutputFile("\n");
	    
    // Line 2 : 
    // write the cable diameter to the Msg parameter
    S << property[i]->Diam.value;
    Msg.WriteDataToOutputFile("    ");
    Msg.WriteDataToOutputFile( property[i]->Diam.name );
    // Msg.WriteDataToOutputFile(":     ");
    Msg.WriteDataToOutputFile("         [m]    : ");
    Msg.WriteDataToOutputFile( S.str() );
    Msg.WriteDataToOutputFile("\n");    
    S.str("");S.clear();

    // Line 3 : 
    // write the cable density to the Msg parameter
    S << property[i]->MassDenInAir.value;
    Msg.WriteDataToOutputFile("    ");
    Msg.WriteDataToOutputFile( property[i]->MassDenInAir.name );
    Msg.WriteDataToOutputFile(" [kg/m] : ");
    Msg.WriteDataToOutputFile( S.str() );
    Msg.WriteDataToOutputFile("\n");    
    S.str("");S.clear();

    // Line 5 : 
    // write the cable Young's modulus to the Msg parameter
    S << property[i]->EA.value;
    Msg.WriteDataToOutputFile("    ");
    Msg.WriteDataToOutputFile( property[i]->EA.name );
    Msg.WriteDataToOutputFile("           [kN]   : ");
    Msg.WriteDataToOutputFile( S.str() );
    Msg.WriteDataToOutputFile("\n");    
    S.str("");S.clear();

    // Line 6 : 
    // write the cable friction coefficient to the Msg parameter
    S << property[i]->CB.value;
    Msg.WriteDataToOutputFile("    ");
    Msg.WriteDataToOutputFile( property[i]->CB.name );
    Msg.WriteDataToOutputFile("                  : ");
    Msg.WriteDataToOutputFile( S.str() );
    Msg.WriteDataToOutputFile("\n\n");    
    S.str("");S.clear();
  }
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
void MAP_OtherStateType_class::
writeEnvironmentData  ( MAP_Message_class &Msg )
{
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
 * WriteNodeData
 *
 * @input : Msg          -- Error message status
 * ============================================================================
 */
void MAP_OtherStateType_class::
WriteNodeData( MAP_Message_class &Msg )
{
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

    if( j+FOUR > node.size( ) ) {
      num = node.size() - j;
    } else {
      num = FOUR;
    };
	
    for ( unsigned int col=j ; col<j+num ; col++) {
      // This creates the header for the Node data portion of the output message
      write.clear();
	    
      S << col+1;
      // write += "   Node ";
      write += "Node ";
      write += S.str();
      write += " Data:";
      S.str("");S.clear();
	    
      // To be sure the length of 'write' is not greater than
      // length we are allow to across, this assert statement 
      // is raised.
      assert( write.size()<_TEXT_COLOR::STR_LEN );
	    
      // Fill spaces in the 'write' string with white spaces to 
      // occupy area in the output string (to make the output
      // print in alignment)
      do { write += " "; } while( write.size()<_TEXT_COLOR::STR_LEN);

      // Add column to all rows
      cell[0].push_back( write );
      write.clear();

      switch ( node[col]->GetNodeType() ) {
      case No_Definition:
        // write += "   No Definition ";
        write += "No Definition ";
        break;
      case Fix : 
        write += "Fix";
        break;
      case Vessel :
        write += "Vessel";
        break;
      case Connect : 
        write += "Connect";
        break;
      }

      // To be sure the length of 'write' is not greater than
      // length we are allow to across, this assert statement 
      // is raised. This only because visiable if the condition
      // is violated
      assert( write.size()<_TEXT_COLOR::STR_LEN );

      // Fill spaces in the 'write' string with white spaces to 
      // occupy area in the output string (to make the output
      // print in alignment)
      do{ write += " "; } while( write.size()<_TEXT_COLOR::STR_LEN );
	
      // Add column to all rows
      cell[1].push_back( write );
      write.clear();
	
      // this block pushes strings to each cell in the output string
      // we are writting to
      cell[2].push_back( VarType::WriteNodeData( node[col]->X ) );
      cell[3].push_back( VarType::WriteNodeData( node[col]->Y ) );
      cell[4].push_back( VarType::WriteNodeData( node[col]->Z ) );
      cell[5].push_back( VarType::WriteNodeData( node[col]->M ) );
      cell[6].push_back( VarType::WriteNodeData( node[col]->B ) );
      cell[7].push_back( VarType::WriteNodeData( node[col]->FX ) );
      cell[8].push_back( VarType::WriteNodeData( node[col]->FY ) );
      cell[9].push_back( VarType::WriteNodeData( node[col]->FZ ) );
    }

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
    write += "FX [kN]  :  |  " ;
    for( unsigned int i=j ; i<j+num ; i++) write += cell[7][i];
    write += "\n";

    // Write FY data to message string    
    write += "FY [kN]  :  |  " ;
    for( unsigned int i=j ; i<j+num ; i++) write += cell[8][i]; 
    write += "\n";
  
    // Write FZ data to message string
    write += "FZ [kN]  :  |  " ;
    for( unsigned int i=j ; i<j+num ; i++) write += cell[9][i]; 
    write += "\n\n\n";

    // Finalize output by writting the string to the Msg object
    Msg.WriteDataToOutputFile( write );
  }
};


// ============================================================================
// addNode  >>  Overloaded Member Function
//
// @input : T[0]         --  name of variable
// @input : T[1]         --  X
// @input : T[2]         --  Y
// @input : T[3]         --  Z
// @input : T[4]         --  M
// @input : T[5]         --  B
// @input : T[6]         --  FX
// @input : T[8]         --  FY
// @input : T[9]         --  FZ
// @input : index        -- 
// @input : Msg          -- Error message status
// @input : Error        -- Error code
// ============================================================================
void MAP_OtherStateType_class::
addNode(const std::vector<std::string> &T     ,
        const int                      index  ,
        MAP_ErrStat_class                    &Error ,
        MAP_Message_class                    &Msg   ) 
{  
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
  node_ptr->SetVarType( T[2] , "X"  , index , &Node::X , Error, Msg );
  node_ptr->SetVarType( T[3] , "Y"  , index , &Node::Y , Error, Msg );

  // the case for "Z" is handled a few lines below
  node_ptr->SetVarType( T[5] , "M"  , index , &Node::M , Error, Msg );
  node_ptr->SetVarType( T[6] , "B"  , index , &Node::B , Error, Msg );
  node_ptr->SetVarType( T[7] , "FX" , index , &Node::FX, Error, Msg );
  node_ptr->SetVarType( T[8] , "FY" , index , &Node::FY, Error, Msg );
  node_ptr->SetVarType( T[9] , "FZ" , index , &Node::FZ, Error, Msg );

  // This if-statement evaluates whether the Z variable in the MAP input file reads 
  // as 'depth'. If it does, then we use the depth marked in the calling program. 
  if( boost::iequals( T[4] , "DEPTH" ) ){
    // set the Z node displacement to the default water depth
    node_ptr->SetVarType( boost::lexical_cast<std::string>( this->GetDepth() ) , 
                          "Z"                                                  , 
                          index                                                , 
                          &Node::Z                                             , 
                          Error                                                , 
                          Msg                                                  );
  } else {
    node_ptr->SetVarType( T[4] , "Z"  , index , &Node::Z , Error, Msg );
  };// END if

  // Now psuh the node we just created to the stack
  node.push_back( node_ptr );	

  // Check if Z is smaller than the water depth. This cannot happen. Warn the users they
  // used a potentially wrong value in the MAP input file.
  if ( node_ptr->GetVarTypeValue( &Node::Z ) < this->GetDepth() ){
    node_ptr->SetVarTypeValue( &Node::Z , this->GetDepth() );    
    throw MAP_ERROR_16; // issue error
  };// END if
};


// ============================================================================
// addElement  
// ============================================================================
void MAP_OtherStateType_class::
addElement( const std::vector<std::string> &T     ,
            const int                      index  ,
            MAP_ErrStat_class                    &Error ,
            MAP_Message_class                    &Msg   ) 
{
  // create new memory on the heap (dynamically allocated)
  Element_ptr element_ptr( new Element );
  std::string error_output = "";

  element_ptr->SetLineProperty( T[1] , property, Error, Msg);
  
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
    element_ptr->SetFairlead( top , node );

    // convert string to integer
    bottom = boost::lexical_cast<int>( T[3] );

    // make sure 'top' is not out of index
    if ( bottom > node.size( ) ) throw MAP_ERROR_22;

    // 'node' is a vector of nodes. 'bottom' is
    // the index we are setting
    element_ptr->SetAnchor  ( bottom , node );
    
    // set element unstretched length
    // Here we set the Lu variable in class Element
    element_ptr->SetLu( T[2], index, Error, Msg );

    // Finally, set H and V as iterated value
    element_ptr->SetH_and_V_flags(); 
    element_ptr->SetH( index , Error , Msg ); // this calls VarType::SetGenericVarType
    element_ptr->SetV( index , Error , Msg ); // this calls VarType::SetGenericVarType

    element_ptr->SetHX(); // purely for writing Hx data to the output file
    element_ptr->SetHY(); // purely for writing Hy data to the output file

    int cnt = 0;
    if( element_ptr->GetLuFlag()==false ) cnt++;
    if( element_ptr->GetHFlag() ==false ) cnt++;
    if( element_ptr->GetVFlag() ==false ) cnt++;

    // if too many variables are iterated, issue a warning
    if (cnt >=3 ) throw MAP_ERROR_24;

    // if too little variables are iterated, issue a warning
    if (cnt <=1 ) throw MAP_ERROR_25;

    for (unsigned int i=5 ; i<T.size()-1 ; i++ ) {
      this->setElementOptions( element_ptr , T[i] , Error, Msg );
    }

    // This last step creates the node. This is a boost::shared_ptr
    element.push_back( element_ptr );
  } catch ( boost::bad_lexical_cast const& ) {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_23 );
    str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_23 );

    Error.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordToErrorList( str );

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
    Msg.RecordToErrorList( str );

    // create context to write error in summary.map        
    str.erase(2,1);                                        
    S.str("");S.clear();                                   
    S << ">>>> " +str +"    :    MAP error in file ";      
    S << __FILE__;                                         
    S << " on line ";                                      
    S << __LINE__;                                         
    S << " occurred.";                                     
    Msg.WriteErrorToOutputFile( S.str() );                            
  }
};


// ============================================================================
// associateElementVarTypeToNWTCType
//
// @input : i            -- 
// @input : P            -- 
// @input : C            --   
//
// @todo : V and H should be local states instead of assigned to parameters
//         (since they change in time)
// ============================================================================
void MAP_OtherStateType_class::
associateElementVarTypeToNWTCType( const int                     index ,
                                   MAP_ParameterType_class       &P    , 
                                   MAP_ConstraintStateType_class &C    ,
                                   MAP_ErrStat_class                   &Err  ,
                                   MAP_Message_class                   &Msg  ) 
{
  // Assign element H variable
  if ( this->elementVarTypeBool( index , &Element::H )==true ) {
    this->setElementReferenceToMAPType( index , P , &Element::H );  
  } else {                                                                    
    this->setElementReferenceToMAPType( index , C , &Element::H );  
  }
                                                                               	
  // Assign element V variable
  if ( this->elementVarTypeBool( index , &Element::V )==true ) {              
    this->setElementReferenceToMAPType( index , P , &Element::V );     
  } else {                                                                    
    this->setElementReferenceToMAPType( index , C , &Element::V );  
  }
                                                                              
  // Assign element Lu variable
  if ( this->elementVarTypeBool( index , &Element::Lu )==true ){ 
    this->setElementReferenceToMAPType( index , P , &Element::Lu );    
  } else {                                                                    
    this->setElementReferenceToMAPType( index , C , &Element::Lu ); 
  }
};


// ============================================================================
// checkElementVarTypeReferences
//
// As a last step, we check the reference_counter for all VarType (in Node 
// Element and CableLibrary) to ensure it is pointing to one MAP_ParameterType_class,               
// MAP_ConstraintStateType_class or MAP_InputType_class.
//
// @input : i            -- 
// @input : Msg          -- Error message status
// @input : Error        -- Error code
// ============================================================================
void MAP_OtherStateType_class::
checkElementVarTypeReferences( const int   index  , 
                               MAP_ErrStat_class &Error , 
                               MAP_Message_class &Msg   )
{
  if ( this->checkElementReference( index , &Element::Lu ) != 1 ){    
    std::string str = "";
    std::ostringstream S;
    S << index+1;
    str += boost::lexical_cast < std::string > ( MAP_ERROR_42 );
    str += "] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_42 ) + S.str() + ".";

    Error.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordToErrorList( str );
  }
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
void MAP_OtherStateType_class::
checkNodeVarTypeReferences( const int index    , 
                            MAP_ErrStat_class &Error , 
                            MAP_Message_class &Msg   )
{
  try {
    if ( this->checkNodeReference( index , &Node::X ) != 1 ){
      throw MAP_ERROR_43;
    }

    if ( this->checkNodeReference( index , &Node::Y ) != 1 ){
      throw MAP_ERROR_44;
    }

    if ( this->checkNodeReference( index , &Node::Z ) != 1 ){    
      throw MAP_ERROR_45;
    }

    if ( this->checkNodeReference( index , &Node::M ) != 1 ){    
      throw MAP_ERROR_46;
    }

    if ( this->checkNodeReference( index , &Node::B ) != 1 ){    
      throw MAP_ERROR_47;
    }

    if ( this->checkNodeReference( index , &Node::FX ) != 1 ){   
      throw MAP_ERROR_48;
    }

    if ( this->checkNodeReference( index , &Node::FY ) != 1 ){    
      throw MAP_ERROR_49;
    }

    if ( this->checkNodeReference( index , &Node::FZ ) != 1 ){    
      throw MAP_ERROR_50;
    }
  } catch( MAP_ERROR_CODE &code ) {
    std::string str = "";
    std::ostringstream S;
    S << index+1;
    str += boost::lexical_cast < std::string > ( code );
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";
    Error.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordToErrorList( str );
                                                          
    // create context to write error in summary.map  
    str.erase(2,1);                                 
    S.str("");S.clear();                            
    S << ">>>> " +str +"    :    MAP error in file "; 
    S << __FILE__;                                   
    S << " on line ";                               
    S << __LINE__;                                  
    S << " occurred.";                              
    Msg.WriteErrorToOutputFile( S.str() );                     
  }
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
void MAP_OtherStateType_class::
setElementOptions( Element_ptr       &P     , 
                   const std::string &T     , 
                   MAP_ErrStat_class       &Error , 
                   MAP_Message_class       &Msg   )
{
  EnumParser <ElementOptions> parser;
  std::string WORD = boost::to_upper_copy( T );

  switch( parser.parseElementOptions( WORD ) ) {
  case PLOT :  // plot cable using matplotlib
    P->SetPlotFlag( true );
    break;
  case X_POS : // print node X position to output file
    P->SetXPosFlag( true );
    break;
  case Y_POS : // print node Y position to output file
    P->SetYPosFlag( true );
    break;
  case Z_POS : // print node Z position to output file
    P->SetZPosFlag( true );
    break;
  case X_FORCE : // print node sum forces in the X direction to output file
    P->SetXForceFlag( true );
    break;
  case Y_FORCE : // print node sum forces in the Y direction to output file
    P->SetyForceFlag( true );
    break;
  case Z_FORCE : // print node sum forces in the Z direction to output file
    P->SetZForceFlag( true );
    break;
  case LINE_TENSION : // print tension along the line to output file
    P->SetLineTensionFlag( true );
    break;
  case OMIT_CONTACT : // ignore cable/seabed contact
    P->SetOmitContactFlag( true );
    break;
  case LAY_LENGTH : // length of line touching the seafloor
    P->SetLayLengthFlag( true );
    break;
  case FAIR_TENSION :
    P->SetFairTensionFlag( true );
    break;
  case ANCH_TENSION :
    P->SetAnchTensionFlag( true );
    break;
  default :
    // If the string is an empy space ' ', then avoid printing an error to the screen
    if( boost::iequals( T , " " ) == false ){
      Error.set_error_key( MAP_WARNING );
      Msg.RecordToWarningList("Could not interpret '" 
                              + T 
                              + "' as a valid element option in the MAP input file.");
    }
    break;
  }
};


// ============================================================================
// writeElementData
//
// @input : Msg          -- Error message status
// ============================================================================
void MAP_OtherStateType_class::
writeElementData( MAP_Message_class &Msg )
{
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
                               + element[i]->GetElementName( )
                               + "\n");
	
    // Line 2: Write X position for the element
    // set the first and last node position 
    first << element[i]->GetAnchorPosition( &Node::X );
    last << element[i]->GetFairleadPostion( &Node::X );

    // this write the data to a message status string
    writeXYZData( "    X Position [m] : |  "       , 
                  first.str()                       , 
                  element[i]->GetAnchorFlag( &Node::X ),
                  last.str()                        , 
                  element[i]->GetFairleadFlag( &Node::X ),
                  Msg );

    // clear ostringstream for next use
    first.str(""); first.clear();
    last.str(""); last.clear();

    // Line 3: Write Y position for the element
    // set the first and last node position 
    first << element[i]->GetAnchorPosition( &Node::Y );
    last << element[i]->GetFairleadPostion( &Node::Y );

    // this write the data to a message status string
    writeXYZData( "    Y Position [m] : |  "       , 
                  first.str()                       , 
                  element[i]->GetAnchorFlag( &Node::Y ),
                  last.str()                        , 
                  element[i]->GetFairleadFlag( &Node::Y ),
                  Msg );

    // clear ostringstream for next use
    first.str(""); first.clear();
    last.str(""); last.clear();

    // Line 4: Write Z position for the element
    // set the first and last node position 
    first << element[i]->GetAnchorPosition( &Node::Z );
    last << element[i]->GetFairleadPostion( &Node::Z );

    // this write the data to a message status string
    writeXYZData( "    Z Position [m] : |  "       , 
                  first.str()                       , 
                  element[i]->GetAnchorFlag( &Node::Z ),
                  last.str()                        , 
                  element[i]->GetFairleadFlag( &Node::Z ),
                  Msg );

    // clear ostringstream for next use
    first.str(""); first.clear();
    last.str(""); last.clear();

    // Line 5: write the unstretched element length
    S << element[i]->GetLu();
    if ( element[i]->GetLuFlag() == false ) write +="(";// _TEXT_COLOR::BLUE;
    write += S.str();
	
    // Remove the blue text option
    if ( element[i]->GetLuFlag() == false ) write +=")";//write += _TEXT_COLOR::END;

    Msg.WriteDataToOutputFile( "    Length     [m] : |  "
                        + write
                        + "\n");
    write.clear();
    S.str(""); S.clear();


    // Line 6: write the H horizontal element force
    S << element[i]->GetH();
    if ( element[i]->GetHFlag() == false ) write += "(";//_TEXT_COLOR::BLUE;
    write += S.str();
	
    // Remove the blue text option
    if ( element[i]->GetHFlag() == false ) write += ")";//_TEXT_COLOR::END;

    Msg.WriteDataToOutputFile( "    H          [kN]: |  "
                        + write
                        + "\n");
    write.clear();
    S.str(""); S.clear();	

    // Line 6: write V vertical element force
    S << element[i]->GetV();
    if ( element[i]->GetVFlag() == false ) write += "(";//_TEXT_COLOR::BLUE;
    write += S.str();
	
    // Remove the blue text option
    if ( element[i]->GetVFlag() == false ) write += ")";//write += _TEXT_COLOR::END;

    Msg.WriteDataToOutputFile( "    V          [kN]: |  "
                        + write
                        + "\n\n\n");
    write.clear();
    S.str(""); S.clear();	
  }
};


// ============================================================================
// writeXYZData
//
// @input : position     -- 
// @input : tail         -- 
// @input : tail_bool    -- 
// @input : head         -- 
// @input : head_bool    -- 
// @input : Msg          -- Error message status
// ============================================================================
void MAP_OtherStateType_class::
writeXYZData( const std::string &position , 
              const std::string &tail     ,
              const bool        tail_bool ,
              const std::string &head     ,
              const bool        head_bool ,
              MAP_Message_class &Msg      )
{    
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


// ============================================================================
// summary
//
// @output : 
// ============================================================================
std::string MAP_OtherStateType_class::
summary( ) 
{
  // Now print the CableLibrary, Node and Element
  // variables to the Msg parameter

  MAP_Message_class T;

  T.WriteDataToOutputFile("\n");

  this->writeEnvironmentData          ( T );
  this->writeCableLibraryData         ( T );
  this->WriteNodeData                 ( T );
  this->writeElementData              ( T );

  if( numeric_method.GetMsqsKFlag() ) { 
    checkpoint();
    // Disable the ability to write the linearized stiffness matrix for now
    // until this is verified to work correctly
    //
    // If the numerics routine is initialized, run the 
    // stiffness matrix linearization
    // @todo: find a better/robust/accurate way of calculating the 
    //        linearized stiffness matrix.
    if ( this->isNumericsUninitialized()==false ){ 
      this->writeLinearizedStiffnessMatrix( T );
    }
  }

  return T.GetDataString();
};


// ============================================================================
// setCableLibraryReference
//   
// This function is always called to assign CableLibrary as input in the 
// MAP_ParameterType_class object
//
// @input : T            -- 
// @input : i            -- 
// ============================================================================
void MAP_OtherStateType_class::
setCableLibraryReference( MAP_ParameterType_class &T    , 
                          const int               index ) 
{ 
  // This pushes a reference of all CableLibrary variables into MAP_ParameterType_class
  T.pushVar( &property[index]->Diam );         // cable diameter
  T.pushVar( &property[index]->MassDenInAir ); // cable mass density
  T.pushVar( &property[index]->EA );           // cable stiffness
  T.pushVar( &property[index]->CB );           // cable/seabed friction coefficient
};


// ============================================================================
// checkNodeReference
//
// @input : i            -- 
// @input : ptr          -- 
//
// @output : 
// ============================================================================
int MAP_OtherStateType_class::
checkNodeReference( const int       index , 
                    VarType Node::* ptr   ) const 
{ 
  // @todo : remove assertion in next version of MAP
  assert( ((*node[index]).*ptr).reference_counter==1  ) ;
  return ((*node[index]).*ptr).reference_counter; 
};


// ============================================================================
// getList
//
// This plots the "details" of the MAP_Other type in Python
// ============================================================================
#ifdef WITH_PYTHON
std::string MAP_OtherStateType_class::
getList( )
{
  return this->OtherStateDataTypes.list();
}
#endif

// ============================================================================
// checkElementReference
//
// @input : i            -- 
// @input : prt          -- 
//
// @output :
// @todo : remove assertion in next version of MAP (should write error message) 
// ============================================================================
int MAP_OtherStateType_class::
checkElementReference( const int          index , 
                       VarType Element::* ptr   ) const 
{ 
  assert( ((*element[index]).*ptr).reference_counter==1  ) ;
  return ((*element[index]).*ptr).reference_counter; 
};


// ============================================================================
// addDepth
//
// @input : T            -- 
// @input : Msg          -- Error message status
// @input : Error        -- Error code
// ============================================================================
void MAP_OtherStateType_class::
addDepth( const std::string &T     ,
          MAP_ErrStat_class       &Error ,
          MAP_Message_class       &Msg   ) 
{    
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
    Msg.RecordToErrorList( str );
                                                                                                     
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
  }
};


// ============================================================================
// addGravity 
//
// @input : T            -- 
// @input : Msg          -- Error message status
// @input : Error        -- Error code
// ============================================================================
void MAP_OtherStateType_class::
addGravity( const std::string &T     ,
            MAP_ErrStat_class       &Error ,
            MAP_Message_class       &Msg   ) 
{    
  if ( T == "" ) throw MAP_ERROR_8;

  // make sure the variable can be converted to a double
  try {
    this->gravity = boost::lexical_cast<double>( T );
  } catch ( boost::bad_lexical_cast const& ) {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_11 );
    str += " ] : " + MAP_ERROR_CODE_to_string.at( MAP_ERROR_11 );
        
    Error.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordToErrorList( str );
                                                        
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
  }
};


// ============================================================================
// addSeaDensity
//
// @input : T            -- 
// @input : Msg          -- Error message status
// @input : Error        -- Error code
// ============================================================================
void MAP_OtherStateType_class::
addSeaDensity( const std::string &T     ,
               MAP_ErrStat_class       &Error ,
               MAP_Message_class       &Msg   ) 
{  
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
    Msg.RecordToErrorList( str );

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
  }
};


// ============================================================================
// setNodeReferenceToUserData
//
// Passes the address of a 'Connect' node to UserData so the Newton force 
// balance equation can be access in the PETSc solving function
//
// @todo : a case may come up where one component of Newton's equations does not have
//         to be solved.  When this comes us, we include a better method to count 
//         the number of Newton equations solved.
//
// @input : i            -- 
// ============================================================================
void MAP_OtherStateType_class::
setNodeReferenceToUserData( const int index )
{ 
  Node* raw_Node = node[index].get();
  this->user_data.pushNodeBack( raw_Node ); 
};


// ============================================================================
// setMessageReferenceToUserData
//
// @input :              -- 
// ============================================================================
void MAP_OtherStateType_class::
setMessageReferenceToUserData( MAP_Message_class &Msg ) 
{ 
  this->user_data.setMessage( Msg );
};


// ============================================================================
// setErrorStatusReferenceToUserData
//
// @input :              -- 
// ============================================================================
void MAP_OtherStateType_class::
setErrorStatusReferenceToUserData( MAP_ErrStat_class &error ) 
{ 
  this->user_data.setErrorCode( error );
};


// ============================================================================
// setElementReferenceToUserData
//
// Passes the adress of each element so that the continuous analytical catenary
// equation can be accessed in the PETSc solve functions. 
//
// @todo : a case may come up where we may need only one (instead of two)
//         catenary equations to solve the problem. This will come up when H, Lu
//         is defined, but not V. In the future, write logic into the program so
//         that the appropriate number of catenary equations are used to solve
//         the problem. 
//
// @input : i            -- 
// ============================================================================
void MAP_OtherStateType_class::
setElementReferenceToUserData( int i )
{ 
  Element* raw_Element = element[i].get();
  this->user_data.pushElementBack( raw_Element );  
};


// ============================================================================
// setMAP_ConstraintStateType_classReferenceToUserData
//
// Passes the address of the constraint states (the variables we are iterating/
// solving)
//
// @input : T            -- 
// ============================================================================
void MAP_OtherStateType_class::
setMAP_ConstraintStateType_classReferenceToUserData( MAP_ConstraintStateType_class &T )
{ 
  this->user_data.initializeConstriant( T );
};


// ============================================================================
// plot
//
// User option to plot the 3D cable profile with MatPlotLib
//
// @input : Error        -- Error code
// @input : Msg          -- Error message status
// ============================================================================
#ifdef WITH_PYTHON
void MAP_OtherStateType_class::
plot( MAP_ErrStat_class &Error , 
      MAP_Message_class &Msg   ) 
{ 
  Msg.MessageClean();

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
      flag  = element[i]->GetElementPlotArray(x,y,z,Error,Msg);

      if ( flag == true ) {      
        // get color
        if ( element[i]->GetElementName( ) == "nylon" ) {
          //color = ",color=[.725,.725,.113],lw=3";  // dark yello
          color = ",color=[0,.4,.8],lw=3"; // light blue
          //color = ",color='b',lw=3";
          //color = ",color=[.784,.568,.568],lw=3"; // pink
          //color = ",color=[.380,.709,.811],lw=3"; // blue
          //color = ",color=[0.015,0.749,1.0],lw=3"; // blue
          //color = ",color=[.898,.537,0],lw=3"; // orange
        } else if ( element[i]->GetElementName( ) == "nylon2" ) {
          //color = ",color=[.282,.487,.487],lw=3";
          color = ",color=[.917,.917,.674],lw=3"; // yello
          //color = ",color=[0,.7,1],lw=3";
        } else if (element[i]->GetElementName( ) == "steel" ) {
          //color = ",color=[.725,.725,.113],lw=3";  // dark yello
          color = ",color=[1,102/255,102/255],lw=3";//color = ",color='b',lw=3";
          //color = ",color=[.725,.725,.113],lw=3";  // dark yello
          //color = ",color='r',lw=3";
        } else if (element[i]->GetElementName( ) == "chainup" ) {
          color = ",color='k',lw=3";
        } else if (element[i]->GetElementName( ) == "chain" ) {
          color = ",color='k',lw=3";
          //color = ",color=[0.015,0.749,1.0],lw=3"; // blue
        } else if (element[i]->GetElementName( ) == "polyester" ) {
          color = ",color=[0.0,0.305,0.478],lw=3";
          //color = ",color=[.725,.725,.113],lw=3";  // dark yello
        } else { 
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
        out += color;
        out += ");";

        // =======  Plot black cirecly on line  ======       <--------------------------+
        std::vector <std::string> temp_string_x;             
        std::vector <std::string> temp_string_y;             
        std::vector <std::string> temp_string_z;             
        
        boost::split(temp_string_x,x,boost::is_any_of(" , "));
        boost::split(temp_string_y,y,boost::is_any_of(" , "));
        boost::split(temp_string_z,z,boost::is_any_of(" , "));
        
        out += "ax.plot([";						
        out += temp_string_x[0] + " , " + temp_string_x[300] + "],[";
        out += temp_string_y[0] + " , " + temp_string_y[300] + "],[";
        out += temp_string_z[0] + " , " + temp_string_z[300] + "],";
        out += "'ko')";						
        //============== <END> =========================================================

        plot_string.push_back( out );
	
      } else {
        std::string str = "";
        MAPSetUniversalErrorStat( MAP_ERROR_81 , str, Error, Msg );
      }
    }
  }
    
  if ( plot_string.size() != 0 ) {    
    Py_Initialize();
    PyRun_SimpleString("import matplotlib as mpl");
    PyRun_SimpleString("import math");
    
//PyRun_SimpleString("from matplotlib import rc");
//PyRun_SimpleString("rc('text',usetex=True)");
//PyRun_SimpleString("rc('font',family='serif')");

    //PyRun_SimpleString("rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})");
    //PyRun_SimpleString("rc('text', usetex=True)");
    PyRun_SimpleString("from mpl_toolkits.mplot3d import Axes3D");
    PyRun_SimpleString("import numpy as np");
    PyRun_SimpleString("import matplotlib.pyplot as plt");
    PyRun_SimpleString("fig = plt.figure()");
    PyRun_SimpleString("ax = Axes3D(fig)");

//// Make background black
//PyRun_SimpleString("ax.set_axis_bgcolor([0.13725,0.13725,0.13725])");

    // Now include the line points from Element::GetElementPlotArray( ... )
    for ( unsigned int i=0 ; i<plot_string.size() ; i++ ){
      PyRun_SimpleString( plot_string[i].c_str() );
    }
    
//    // Give axis limits
    //PyRun_SimpleString("ax.set_xlim3d(0, 400)");
    //PyRun_SimpleString("ax.set_ylim3d(-40, 40)");
    //PyRun_SimpleString("ax.set_zlim3d(-350, 100)");

//// Make axis tick white
//PyRun_SimpleString("for t in ax.w_xaxis.get_ticklines(): t.set_color('w')");
//PyRun_SimpleString("for t in ax.w_yaxis.get_ticklines(): t.set_color('w')");
//PyRun_SimpleString("for t in ax.w_zaxis.get_ticklines(): t.set_color('w')");
//PyRun_SimpleString("for tick in ax.w_xaxis.get_major_ticks() + ax.w_yaxis.get_major_ticks() + ax.w_zaxis.get_major_ticks():\n for child in tick.get_children(): \n     child.set_color('w')");
//PyRun_SimpleString("ax.set_xlabel(r'X [m]', color='w',fontsize=16)");
//PyRun_SimpleString("ax.set_ylabel(r'Y [m]', color='w',fontsize=16)");
//PyRun_SimpleString("ax.set_zlabel(r'Depth [m]', color='w',fontsize=16)");
//PyRun_SimpleString("myAXINFO = { 'x': {'i': 0, 'tickdir': 1, 'juggled': (1, 0, 2), 'color': (1.0, 1.0, 1.0, 1.0, 1.0)}, 'y': {'i': 1, 'tickdir': 0, 'juggled': (0, 1, 2), 'color': (1.0, 1.0, 1.0, 1.0, 1.0)}, 'z': {'i': 2, 'tickdir': 0, 'juggled': (0, 2, 1), 'color': (1.0, 1.0, 1.0, 1.0, 1.0)} } ");
//PyRun_SimpleString("ax.w_zaxis._AXINFO = myAXINFO ");
//PyRun_SimpleString("ax.w_yaxis._AXINFO = myAXINFO ");
//PyRun_SimpleString("ax.w_xaxis._AXINFO = myAXINFO ");

    //PyRun_SimpleString("ax.set_xlabel(r'X Displacement [m]', color='k',fontsize=14)");
    //PyRun_SimpleString("ax.set_ylabel(r'Y Displacement [m]', color='k',fontsize=14)");
    //PyRun_SimpleString("ax.set_zlabel(r'Depth Below MSL [m]', color='k',fontsize=14)");
    //PyRun_SimpleString("ax.set_xlabel('X Displacement [m]', color='k',fontsize=14)");
    //PyRun_SimpleString("ax.set_ylabel('Y Displacement [m]', color='k',fontsize=14)");
    //PyRun_SimpleString("ax.set_zlabel('Depth Below MSL [m]', color='k',fontsize=14)");

//// Plotting preferences and save the file
////PyRun_SimpleString("ax.view_init(90, 0)");
//PyRun_SimpleString("ax.view_init(30, 120)");
//PyRun_SimpleString("ax.patch.set_alpha(0.0)");
//PyRun_SimpleString("fig.savefig('/media/sf_Resume/2H/Presentation/Figures/bridle2.pdf', transparent=True)");
	
    // setup axis so a minimum limit is set. 
    PyRun_SimpleString("limit_x = ax.get_xlim3d()");
    PyRun_SimpleString("if math.fabs( (limit_x[0]-limit_x[1])<1 ) : ax.set_xlim3d( (limit_x[0]-1) , (limit_x[1]+1) )");
    PyRun_SimpleString("limit_y = ax.get_ylim3d()");
    PyRun_SimpleString("if math.fabs( (limit_y[0]-limit_y[1])<1 ) : ax.set_ylim3d( (limit_y[0]-1) , (limit_y[1]+1) )");
    PyRun_SimpleString("limit_z = ax.get_zlim3d()");
    PyRun_SimpleString("if math.fabs( (limit_z[0]-limit_z[1])<1 ) : ax.set_zlim3d( (limit_z[0]-1) , (limit_z[1]+1) )");

    // show plot
    PyRun_SimpleString("plt.show()");

    //Py_Exit(0);
  }
};
#endif

void MAP_OtherStateType_class::
GetPyPlotArray( int                 &index ,
                std::vector<double> &X     , 
                std::vector<double> &Y     , 
                std::vector<double> &Z     , 
                MAP_ErrStat_class   &Error , 
                MAP_Message_class   &Msg   ) 
{
  element[index]->GetElementPlotArray(X,Y,Z,Error,Msg);
}


// ============================================================================
// plot
//
// User option to plot the 3D cable profile with MatPlotLib
//
// @input : Error        -- Error code
// @input : Msg          -- Error message status
// ============================================================================
std::vector <std::string> MAP_OtherStateType_class::
plotString( MAP_ErrStat_class &Error , 
            MAP_Message_class &Msg   ) 
{
  MAP_Message_class T;

  std::vector <std::string> plot_string;

  std::string x     = "";
  std::string y     = "";
  std::string z     = "";
  std::string out   = "";
  bool flag         = true;

  for( int i=0 ; i<this->getSizeOfElement() ; i++ ) {
    // if the element has the option flag PLOT, then plot this element to the screen
    if ( this->getElementOptionFlag( i, &Element::PLOT_flag ) ) {	    
      x     = "";
      y     = "";
      z     = "";
      out   = "";
      flag  = element[i]->GetElementPlotArray(x,y,z,Error,Msg);

      if ( flag != true ){
        std::string str = "";
        MAPSetUniversalErrorStat( MAP_ERROR_81 , str, Error, Msg );
      } else { // add the string arrays into the Matplotlib plot command        
        out += "[" + x;
        plot_string.push_back( x );
        
        out += "],[" + y;
        plot_string.push_back( y );
        
        out += "],[" + z + "]";
        plot_string.push_back( z );
      }
    }
  }

  if ( plot_string.size() != 0 ){    
    // Now include the line points from Element::GetElementPlotArray( ... )
    for ( unsigned int i=0 ; i<plot_string.size() ; i++ ){
      T.WriteDataToOutputFile( plot_string[i].c_str() );
    }
  }
  
  return plot_string;
};

// ============================================================================
// Recieves the string inputs arguments from the 'SOLVER OPTIONS' 
// section of the MAP input file and set the PETSc run-time options 
// accordiningly. This is called repeated in a for-loop proportional
// to the number of options in the MAP input file. 
//
// @param   inputStr  Options input string (1 word only)
// @param   msg       error message status
// @return  int       error code. 0 = no error, 1 = no convergence. 
// ============================================================================
void MAP_OtherStateType_class::
SetSolverOptions( const std::string &inputStr , 
                  MAP_Message_class &msg      )
{
  this->numeric_method.setNumericsOptionsString( inputStr );
//this->numeric_method->setNumericsOptionsString( T );
};


// ============================================================================
// Solve
//
// @input : Msg          -- Error message status
// @input : Error        -- Error code
// ============================================================================
int MAP_OtherStateType_class::
Solve( MAP_ErrStat_class &err , 
       MAP_Message_class &msg )
{  
  int errcheck = 0;
  int i = 0;

  try{
    for( i=0 ; i<this->getSizeOfElement() ; i++ ){      
      element[i]->CheckIfCableIsVertical( );
      element[i]->CheckMaximumLineLength( ); // check if the element is double backing or not      
    }

    // Now solve the problem
    errcheck = numeric_method.PetscSolve( *this, err, msg );
    //numeric_method->PetscSolve( *this, Error, Msg );    
    
    for( i=0 ; i<getSizeOfElement() ; i++ ) {
      element[i]->ResetGhostProperties();
    }
    
    if( errcheck==1 ) throw MAP_ERROR_88;
    
  } catch( MAP_ERROR_CODE &code ) { // catch the double backing problem
    std::string str = "";
    std::ostringstream S;
    S << i+1;
    str = S.str();
    MAPSetUniversalErrorStat( code , str, err, msg );    
    return 1;
  }

  return 0;
};


// ============================================================================
// CheckConvergence
//
// Check the residuals to make sure the solution is actually reached. 
//
// @todo : It might be useful to have the residuals printed to the output
//         file for closer inspection post-simulation
// ============================================================================
int MAP_OtherStateType_class::
CheckResidualConvergence( MAP_ErrStat_class        &err , 
                          MAP_Message_class        &msg )
{
  double tol = numeric_method.GetMSQSTol();
  //double tol = numeric_method->GetMSQSTol();
  
  for ( unsigned int i=0 ; i<this->node.size() ; i++ ) {	
    // solve X direction Newton equation
    if (this->node[i]->GetXNewtonEquationFlag()==true){      
      //std::cout << "x : " << this->node[i]->f_x() << std::endl;
      if (this->node[i]->f_x() >= tol) return 1;//throw MAP_ERROR_69;      
    }
	
    // solve Y direction Newton equation
    if (this->node[i]->GetYNewtonEquationFlag()==true){
      //std::cout << "y : " << this->node[i]->f_y() << std::endl;
      if (this->node[i]->f_y() >= tol) return 1;//throw MAP_ERROR_69;
    }

    // solve Z direction Newton equation
    if (this->node[i]->GetZNewtonEquationFlag()==true){
      //std::cout << "z : " << this->node[i]->f_z() << std::endl;
      if (this->node[i]->f_z() >= tol) return 1;//throw MAP_ERROR_69;
    }
  }

  // For each element, the X,Y,Z catenary equation is solved.  
  for ( unsigned int i=0 ; i<this->element.size() ; i++ ) {
    //std::cout << "H : " << this->element[i]->f_h() << std::endl;
    //std::cout << "V : " << this->element[i]->f_v() << std::endl;
    if (this->element[i]->f_h() >= tol) return 1;//throw MAP_ERROR_69;
    if (this->element[i]->f_v() >= tol) return 1;//throw MAP_ERROR_69;
  }

  //std::cout << "\n\n";
  return 0; // no error
};


// ============================================================================
// initializeNumericSolver
//
// @input : Init         -- 
// @input : Msg          -- Error message status
// ============================================================================
void MAP_OtherStateType_class::
initializeNumericSolver( MAP_InitInputType_class &Init  , 
                         MAP_ErrStat_class             &Error , 
                         MAP_Message_class             &Msg ) 
{
  numeric_method.InitializeSolver( *this , Init , Error, Msg );
  //numeric_method->InitializeSolver( *this , Init , Error, Msg );
};


// ============================================================================
// cleanNumericSolver
//
// Terminates the PETSc solve routines and destroys all PETSc objects.  
//
// @input : Error        -- Error code
// @input : Msg          -- Error message status
// ============================================================================
void MAP_OtherStateType_class::
cleanNumericSolver( MAP_ErrStat_class &Error , 
                    MAP_Message_class &Msg   )
{
  numeric_method.PetscEnd( Error, Msg );
  //numeric_method->PetscEnd( Error, Msg );
};


// ============================================================================
// getNumEquations
//
// @output : 
// ============================================================================
int MAP_OtherStateType_class::
getNumEquations( ) 
{
  return this->num_equations;
};


// ============================================================================
// initializeCableElement
//
// @input : index             -- 
// ============================================================================
void MAP_OtherStateType_class::
initializeCableElement( const int   index , 
                        MAP_ErrStat_class &Error , 
                        MAP_Message_class &Msg ) 
{ 
  this->element[index]->InitializeElement( gravity, rho_sea , Error , Msg );
};


// ============================================================================
// initializeCableNode >
//
// this sets the initial conditions for the Fix and Vessel nodes
//
// @input : index             -- node number 
// ============================================================================
void MAP_OtherStateType_class::
initializeCableNode( const int index )
{
  this->node[index]->InitializeNode( );
};


// ============================================================================
// initializeSumForce
//
// set Node.sum_FX, Node.sum_FY and Node.sum_FZ to zero
// ============================================================================
void MAP_OtherStateType_class::
initializeSumForce( ) 
{
  for( unsigned int i=0 ; i<node.size() ; i++) {
    this->node[i]->SetSumForceToZero();
  }
};


// ============================================================================
// solve_sum_force_equation_in_direction_X >
//
// sets the 'solve_X_Newton_equation' flag in class Node
//
// @input : i            -- 
// @input : T            -- 
// ============================================================================
void MAP_OtherStateType_class::
setSolveSumForceEquationInDirectionX( const int index , 
                                      bool      T     )
{ 
  this->node[index]->SetXNewtonEquationFlag( T ); 
}


// ============================================================================
// solve_sum_force_equation_in_direction_Y
//
// sets the 'solve_Y_Newton_equation' flag in class Node
//
// @input : i            -- 
// @input : T            -- 
// ============================================================================
void MAP_OtherStateType_class::
setSolveSumForceEquationInDirectionY( const int index ,
                                      bool      T     )
{ 
  this->node[index]->SetYNewtonEquationFlag( T ); 
}


// ============================================================================
// solve_sum_force_equation_in_direction_Z
//
// sets the 'solve_Z_Newton_equation' flag in class Node
//
// @input : i            -- 
// @input : T            -- 
// ============================================================================
void MAP_OtherStateType_class::
setSolveSumForceEquationInDirectionZ( const int index , 
                                      bool      T     )
{ 
  this->node[index]->SetZNewtonEquationFlag( T ); 
}


// ============================================================================
// getSumForceEquationFlagInDirectionX
//
// @input : i            -- 
//
// @output : 
// ============================================================================
bool MAP_OtherStateType_class::
getSumForceEquationFlagInDirectionX( const int index )
{ 
  return this->node[index]->GetXNewtonEquationFlag( ); 
}


// ============================================================================
// getSumForceEquationFlagInDirectionY
//
// @input : i            -- 
//
// @output : 
// ============================================================================
bool MAP_OtherStateType_class::
getSumForceEquationFlagInDirectionY( const int index )
{ 
  return this->node[index]->GetYNewtonEquationFlag( ); 
}


// ============================================================================
// getSumForceEquationFlagInDirectionZ
//
// @input : i            -- 
//
// @output : 
// ============================================================================
bool MAP_OtherStateType_class::
getSumForceEquationFlagInDirectionZ( const int index )
{ 
  return this->node[index]->GetZNewtonEquationFlag( ); 
}


// ============================================================================
// GetNumberOfOutputs
//
// Gets the number of output values being written to the output file (for codes
// coupled to MAP, such as FAST)
// ============================================================================
int MAP_OtherStateType_class::
MAPCALL_GetNumberOfOutputs( )
{
  int count = 0;
  
  for ( unsigned int i=0 ; i<element.size() ; i++ ) {
    if ( this->getElementOptionFlag( i , &Element::X_POS_flag)        ) count++;
    if ( this->getElementOptionFlag( i , &Element::Y_POS_flag)        ) count++;
    if ( this->getElementOptionFlag( i , &Element::Z_POS_flag)        ) count++;
    if ( this->getElementOptionFlag( i , &Element::X_FORCE_flag)      ) count++;
    if ( this->getElementOptionFlag( i , &Element::Y_FORCE_flag)      ) count++;
    if ( this->getElementOptionFlag( i , &Element::Z_FORCE_flag)      ) count++;
    if ( this->getElementOptionFlag( i , &Element::LAY_LENGTH_flag)   ) count++;
    if ( this->getElementOptionFlag( i , &Element::LINE_TENSION_flag) ) count+=10;
    if ( this->getElementOptionFlag( i , &Element::FAIR_TENSION_flag) ) count++;
    if ( this->getElementOptionFlag( i , &Element::ANCH_TENSION_flag) ) count++;
  }

  return count;
};



// ============================================================================
// MAPCALL_GetOutputHeaderUnits
//
// @input : i            -- 
// @input : Msg          -- Error message status
// ============================================================================
void MAP_OtherStateType_class::
MAPCALL_GetOutputHeaderUnits( char **arr )
{
  std::string tempStr = "";
  int count = 0;
  int TEN = 10;

  for ( unsigned int i=0 ; i<element.size() ; i++ ) {
    if ( this->getElementOptionFlag( i , &Element::X_POS_flag) ){    
      strcpy( arr[count], "[m]    " );
      count++;      
    }

    if ( this->getElementOptionFlag( i , &Element::Y_POS_flag) ){
      strcpy( arr[count], "[m]    " );
      count++;      
    }

    if ( this->getElementOptionFlag( i , &Element::Z_POS_flag) ){
      strcpy( arr[count], "[m]    " );
      count++;      
    }

    if ( this->getElementOptionFlag( i , &Element::X_FORCE_flag) ){
      strcpy( arr[count], "[kN]   " );
      count++;      
    }

    if ( this->getElementOptionFlag( i , &Element::Y_FORCE_flag) ){
      strcpy( arr[count], "[kN]   " );
      count++;      
    }

    if ( this->getElementOptionFlag( i , &Element::Z_FORCE_flag) ){
      strcpy( arr[count], "[kN]   " );
      count++;      
    }

    if ( this->getElementOptionFlag( i , &Element::LAY_LENGTH_flag) ){
      strcpy( arr[count], "[m]    " );
      count++;      
    }

    if ( this->getElementOptionFlag( i , &Element::LINE_TENSION_flag) ){
      strcpy( arr[count  ], "[kN]   " );
      strcpy( arr[count+1], "[kN]   " );
      strcpy( arr[count+2], "[kN]   " );
      strcpy( arr[count+3], "[kN]   " );
      strcpy( arr[count+4], "[kN]   " );
      strcpy( arr[count+5], "[kN]   " );
      strcpy( arr[count+6], "[kN]   " );
      strcpy( arr[count+7], "[kN]   " );
      strcpy( arr[count+8], "[kN]   " );
      strcpy( arr[count+9], "[kN]   " );
      count+=10;      
    }

    if ( this->getElementOptionFlag( i , &Element::FAIR_TENSION_flag) ){
      strcpy( arr[count], "[kN]   " );
      count++;      
    }

    if ( this->getElementOptionFlag( i , &Element::ANCH_TENSION_flag) ){
      strcpy( arr[count], "[kN]   " );
      count++;      
    }
  }
};


// ============================================================================
// MAPCALL_GetOutputHeader
//
// ============================================================================
void MAP_OtherStateType_class::
MAPCALL_GetOutputHeader( char **arr  ) 
{
  std::string tempStr = "";
  int count = 0;
  int TEN = 10;

  for ( unsigned int i=0 ; i<element.size() ; i++ ) {
    if ( this->getElementOptionFlag( i , &Element::X_POS_flag) ){
      tempStr = VarType::WriteGenericVarType_name( i , element[i]->fairlead->X );
      strcpy( arr[count], tempStr.c_str() );
      count++;
    }
    
    if ( this->getElementOptionFlag( i , &Element::Y_POS_flag) ){
      tempStr = VarType::WriteGenericVarType_name( i , element[i]->fairlead->Y );
      strcpy( arr[count], tempStr.c_str() );
      count++;
    }
    
    if ( this->getElementOptionFlag( i , &Element::Z_POS_flag) ){
      tempStr = VarType::WriteGenericVarType_name( i , element[i]->fairlead->Z );
      strcpy( arr[count], tempStr.c_str() );
      count++;
    }
    
    if ( this->getElementOptionFlag( i , &Element::X_FORCE_flag) ){
      tempStr = VarType::WriteGenericVarType_name( i , element[i]->HX );
      strcpy( arr[count], tempStr.c_str() );
      count++;
    }
    
    if ( this->getElementOptionFlag( i , &Element::Y_FORCE_flag) ){
      tempStr = VarType::WriteGenericVarType_name( i , element[i]->HY );
      strcpy( arr[count], tempStr.c_str() );
      count++;
    }
    
    if ( this->getElementOptionFlag( i , &Element::Z_FORCE_flag) ){
      tempStr = VarType::WriteGenericVarType_name( i , element[i]->V );
      strcpy( arr[count], tempStr.c_str() );
      count++;
    }
    
    if ( this->getElementOptionFlag( i , &Element::LAY_LENGTH_flag) ){
      tempStr = element[i]->GetLayLengthStringHeader(i);
      strcpy( arr[count], tempStr.c_str() );
      count++;
    }
    
    if ( this->getElementOptionFlag( i , &Element::LINE_TENSION_flag) ){
      for( int j=0 ; j<TEN ; j++ ) {
        tempStr = element[i]->GetLineTensionStringHeaderForFast(i,j);        
        strcpy( arr[count], tempStr.c_str() );
        count++;
      }

    }  

    if ( this->getElementOptionFlag( i , &Element::FAIR_TENSION_flag) ){
      tempStr = element[i]->GetFairleadTensionStringHeader(i);
      strcpy( arr[count], tempStr.c_str() );
      count++;
    }

    if ( this->getElementOptionFlag( i , &Element::ANCH_TENSION_flag) ){
      tempStr = element[i]->GetAnchorTensionStringHeader(i);
      strcpy( arr[count], tempStr.c_str() );
      count++;
    }
  }
};


// ============================================================================
// getOutputStreamHeader
//
// @input : i            -- 
// @input : Msg          -- Error message status
// ============================================================================
void MAP_OtherStateType_class::
getOutputStreamHeader( const int   index ,
                       MAP_Message_class &Msg  ) 
{
  // print the simulation time
  if ( is_coupled_to_FAST==false) { // if MAP isn't coupled to FAST, don't print time in the header
    if (index+1 == 1 ){
      this->output_string += "Time           ";
    }
  }

  if ( this->getElementOptionFlag( index , &Element::X_POS_flag) ){
    this->output_string += VarType::WriteGenericVarType_name( index , element[index]->fairlead->X );
  }

  if ( this->getElementOptionFlag( index , &Element::Y_POS_flag) ){
    this->output_string += VarType::WriteGenericVarType_name( index , element[index]->fairlead->Y );
  }

  if ( this->getElementOptionFlag( index , &Element::Z_POS_flag) ){
    this->output_string += VarType::WriteGenericVarType_name( index , element[index]->fairlead->Z );
  }

  if ( this->getElementOptionFlag( index , &Element::X_FORCE_flag) ){
    this->output_string += VarType::WriteGenericVarType_name( index , element[index]->HX );
  }

  if ( this->getElementOptionFlag( index , &Element::Y_FORCE_flag) ){
    this->output_string += VarType::WriteGenericVarType_name( index , element[index]->HY );
  }

  if ( this->getElementOptionFlag( index , &Element::Z_FORCE_flag) ){
    this->output_string += VarType::WriteGenericVarType_name( index , element[index]->V );
  }

  if ( this->getElementOptionFlag( index , &Element::LAY_LENGTH_flag) ){
    this->output_string += element[index]->GetLayLengthStringHeader(index);
  }

  if ( this->getElementOptionFlag( index , &Element::LINE_TENSION_flag) ){
    this->output_string += element[index]->GetLineTensionStringHeader(index);
  }

  if ( this->getElementOptionFlag( index , &Element::FAIR_TENSION_flag) ){
    this->output_string += element[index]->GetFairleadTensionStringHeader(index);;
  }

  if ( this->getElementOptionFlag( index , &Element::ANCH_TENSION_flag) ){
    this->output_string += element[index]->GetAnchorTensionStringHeader(index);;
  }
  
  // end the line
  if ( is_coupled_to_FAST==false) { // if MAP isn't coupled to FAST, don't end the line
    if (index+1 >= this->getSizeOfElement() ){
      this->output_string += "\n";
    }
  }
};

// ============================================================================
//
// ============================================================================
std::string &MAP_OtherStateType_class::
GetOutputString( ) {
  return this->output_string;
}

// ============================================================================
// getOutputStreamUnits
//
// @input : i            -- 
// @input : Msg          -- Error message status
// ============================================================================
void MAP_OtherStateType_class::
getOutputStreamUnits( const int index  , 
                      MAP_Message_class &Msg ) 
{
  // print the simulation time
  if ( is_coupled_to_FAST==false) { // if MAP isn't coupled to FAST, don't print time in the header
    if (index+1 == 1 ){
      this->output_string += "[s]            ";
    }
  } 


  if ( this->getElementOptionFlag( index , &Element::X_POS_flag) ){
    this->output_string += "[m]            " ;
  }

  if ( this->getElementOptionFlag( index , &Element::Y_POS_flag) ){
    this->output_string += "[m]            " ;
  }

  if ( this->getElementOptionFlag( index , &Element::Z_POS_flag) ){
    this->output_string += "[m]            " ;
  }

  if ( this->getElementOptionFlag( index , &Element::X_FORCE_flag) ){
    this->output_string += "[kN]           " ;
  }

  if ( this->getElementOptionFlag( index , &Element::Y_FORCE_flag) ){
    this->output_string += "[kN]           " ;
  }

  if ( this->getElementOptionFlag( index , &Element::Z_FORCE_flag) ){
    this->output_string += "[kN]           " ;
  }

  if ( this->getElementOptionFlag( index , &Element::LAY_LENGTH_flag) ){
    this->output_string += "[m]            " ;
  }

  if ( this->getElementOptionFlag( index , &Element::LINE_TENSION_flag) ){
    this->output_string += "[kN]           " ;
    this->output_string += "[kN]           " ;
    this->output_string += "[kN]           " ;
    this->output_string += "[kN]           " ;
    this->output_string += "[kN]           " ;
    this->output_string += "[kN]           " ;
    this->output_string += "[kN]           " ;
    this->output_string += "[kN]           " ;
    this->output_string += "[kN]           " ;
    this->output_string += "[kN]           " ;
  }

  if ( this->getElementOptionFlag( index , &Element::FAIR_TENSION_flag) ){
    this->output_string += "[kN]           " ;
  }

  if ( this->getElementOptionFlag( index , &Element::ANCH_TENSION_flag) ){
    this->output_string += "[kN]           " ;
  }

  // end the line
  if ( is_coupled_to_FAST==false) { // if MAP isn't coupled to FAST, don't end the line
    if (index+1 >= this->getSizeOfElement() ){
      this->output_string += "\n";
    }
  }
};


// ============================================================================
// getOutputStreamValue
//  
// @input : i            -- 
// @input : Error        -- Error code 
// ============================================================================
void MAP_OtherStateType_class::
getOutputStreamValue( const int   index , 
                      const float time  , 
                      MAP_Message_class &Msg  ) 
{
  // print the simulation time
  if ( is_coupled_to_FAST==false) { // if MAP isn't coupled to FAST, don't end the line
    if (index+1 == 1 ){
      std::ostringstream buff;
      buff << time;
      this->output_string += buff.str();
      do { 
        this->output_string += " ";
      } while( this->output_string.size()<15 ); 
    }
  }

  if ( this->getElementOptionFlag( index , &Element::X_POS_flag) ){
    this->output_string += VarType::WriteGenericVarType_value( index , element[index]->fairlead->X );
  }

  if ( this->getElementOptionFlag( index , &Element::Y_POS_flag) ){
    this->output_string += VarType::WriteGenericVarType_value( index , element[index]->fairlead->Y );
  }

  if ( this->getElementOptionFlag( index , &Element::Z_POS_flag) ){
    this->output_string += VarType::WriteGenericVarType_value( index , element[index]->fairlead->Z );
  }

  if ( this->getElementOptionFlag( index , &Element::X_FORCE_flag) ){
    this->output_string += VarType::WriteGenericVarType_value( index , element[index]->HX );
  }

  if ( this->getElementOptionFlag( index , &Element::Y_FORCE_flag) ){
    this->output_string += VarType::WriteGenericVarType_value( index , element[index]->HY );
  }

  if ( this->getElementOptionFlag( index , &Element::Z_FORCE_flag) ){
    this->output_string += VarType::WriteGenericVarType_value( index , element[index]->V );
  }

  if ( this->getElementOptionFlag( index , &Element::LAY_LENGTH_flag) ){
    this->output_string += element[index]->GetLayLengthString();
  }

  if ( this->getElementOptionFlag( index , &Element::LINE_TENSION_flag) ){
    this->output_string += element[index]->GetLineTensionString();
  }

  if ( this->getElementOptionFlag( index , &Element::FAIR_TENSION_flag) ){
    this->output_string += element[index]->GetFairleadTensionMagnitudeString();
  }

  if ( this->getElementOptionFlag( index , &Element::ANCH_TENSION_flag) ){
    this->output_string += element[index]->GetAnchorTensionMagnitudeString();
  }


  // end the line
  if ( is_coupled_to_FAST==false) { // if MAP isn't coupled to FAST, don't end the line
    if (index+1 >= this->getSizeOfElement() ){
      this->output_string += "\n";
    }
  }
};


// ============================================================================
// openOutputFile
//
// @todo : right now, MAP is hard coded to produce a 'map.out' file.  This is 
//         probably not wanted; it does not allow users to name the output
//         file whatever they want.  
// ============================================================================
void MAP_OtherStateType_class::
openOutputFile( )
{ 
  outputStream->open( "map.out" ); 
};


// ============================================================================
// closeOutputFile
//  
//
// ============================================================================
void MAP_OtherStateType_class::
closeOutputFile( ) 
{ 
  *outputStream << "\n\nMAP successfully terminated";
  outputStream->close();
};


// ============================================================================
// associateNodeVarTypeToNWTCType
//
// @input : i            -- 
// @input : I            -- 
// @input : P            -- 
// @input : C            -- 
// @input : O            -- 
// @input : Msg          -- Error message status
// @input : Error        -- Error code
// ============================================================================
void MAP_OtherStateType_class::
associateNodeVarTypeToNWTCType( const int                     index ,
                                MAP_InputType_class           &I    ,
                                MAP_ParameterType_class       &P    , 
                                MAP_ConstraintStateType_class &C    ,
                                MAP_OutputType_class          &O    ,
                                MAP_ErrStat_class                   &Err  ,
                                MAP_Message_class                   &Msg  )
{  
  // ==========   X Node VarType assignment  ===============     <-------------------------------------------------------------+
  // Assign data in the Node                                                                     
  try {                                                                                                     
    // If X has a fixed value and it IS NOT attached to a vessel,                                           
    // set it as a parameter                                                                                
    if ( this->nodeVarTypeBool( index , &Node::X ) == true && this->GetNodeType( index ) != Vessel ) {      
      // assign as a MAP_ParameterType_class                                                                
      this->setNodeReferenceToMAPType( index , P , &Node::X );                                              
    } else if ( this->nodeVarTypeBool( index , &Node::X ) == true && this->GetNodeType( index ) == Vessel ) { 
      // If X has a fixed value and //it is attached to a vessel, set it as an input assign as a MAP_InputType_class
      this->setNodeReferenceToMAPType( index , I , &Node::X );                                              
    } else if ( this->nodeVarTypeBool( index , &Node::X ) == false ) { // if X is variable (iterated), then it must be a constraint 
      // ensure this node is not a 'Connect' or a 'Vessel'                                                  
      if( this->GetNodeType( index ) == Fix ){                                                              
        // @todo : remove assertion in the next release of MAP                                              
        // once we verify things are working correctly                                                      
        assert( this->GetNodeType( index ) == Connect || this->GetNodeType( index ) == Vessel );            
        throw MAP_ERROR_27;                                                                                 
      }

      // assign as a MAP_ConstraintStateType_class                                                          
      this->setNodeReferenceToMAPType( index , C , &Node::X );                                        

      // Increment the number of equations by                                                         
      // 1 only if the following conditions are                                                       
      // met:                                                                                         
      //   -- X.is_fixed  = false                                                                     
      //   -- FX.is_fixed = true                                                                      
      //   -- Node.type   = Connect                                                                   
      if ( this->nodeVarTypeBool( index, &Node::X )==false  && this->nodeVarTypeBool( index, &Node::FX )==true
           && this->GetNodeType( index ) == Connect ) {                                                       
        // @todo : this is the world's worst way to increment                            
        // the number of equations we are solving. Let's fix                             
        // this before we embarrass ourselves                                            
        this->incrementNumEquations();                                                   
        this->setSolveSumForceEquationInDirectionX( index , true );                      
      } else if ( this->GetNodeType( index ) == Connect ) {                              
        // The node displacement and applied for cannot                                  
        // be simultaneously defined. Here we issue a warning                            
        // message                                                                       
        throw MAP_ERROR_28;                                                              
      }
    } else {                                                                             
      // The node does not meet any initialization criteria                              
      throw MAP_ERROR_29;                                                                
    }
  } catch ( MAP_ERROR_CODE &code ) {                                                     
    std::string str = "";                                                           
    std::ostringstream S;                                                           
    S << index+1;                                                                   
    str += boost::lexical_cast < std::string > ( code );                            
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";            

    Err.set_error_key( MAP_ERROR );                                                 

    // Here is the exception message                                                
    Msg.RecordToErrorList( str );                                                   

    // create context to write error in summary.map                                 
    str.erase(2,1);                                                                 
    S.str("");S.clear();                                                            
    S << ">>>> " +str +"    :    MAP error in file ";                               
    S << __FILE__;                                                                  
    S << " on line ";                                                               
    S << __LINE__;                                                                  
    S << " occurred.";                                                              
    Msg.WriteErrorToOutputFile( S.str() );                                          
  }

  // ==========   Y Node VarType assignment  ===============     <-------------------------------------------------------------+
  // Assign data in the Node                                                                             
  try {                                                                                     
    // If Y has a fixed value and it IS NOT attached to a vessel,                           
    // set it as a parameter                                                                
    if ( this->nodeVarTypeBool( index , &Node::Y ) == true && this->GetNodeType( index ) != Vessel ){    
      // assign as a MAP_ParameterType_class                                                             
      this->setNodeReferenceToMAPType( index , P, &Node::Y );                                            
    } else if ( this->nodeVarTypeBool( index , &Node::Y ) == true && this->GetNodeType( index ) == Vessel ){ 
      // If Y has a fixed value and it is attached to a vessel, set it as an input 

      // assign as a MAP_InputType_class                                                                    
      this->setNodeReferenceToMAPType( index , I , &Node::Y );                                              
    } else if ( this->nodeVarTypeBool( index , &Node::Y ) == false ) { // if Y is variable (iterated), then it must be a constraint 
      // ensure this node is not a 'Connect' or a 'Vessel'                                                  
      if( this->GetNodeType( index ) == Fix ){                                                              
        // @todo : remove assertion in the next release of MAP                                              
        // once we verify things are working correctly                                                      
        assert( this->GetNodeType( index ) == Connect || this->GetNodeType( index ) == Vessel );            
        throw MAP_ERROR_27;                                                                                 
      }

      // assign as a MAP_ConstraintStateType_class                                                          
      this->setNodeReferenceToMAPType( index , C , &Node::Y );                                                

      // Increment the number of equations by                                                                 
      // 1 only if the following conditions are                                                               
      // met:                                                                                                 
      //   -- Y.is_fixed  = false                                                                             
      //   -- FY.is_fixed = true                                                                              
      //   -- Node.type   = Connect                                                                           
      if ( this->nodeVarTypeBool( index, &Node::Y )==false  && this->nodeVarTypeBool( index, &Node::FY )==true
           && this->GetNodeType( index ) == Connect ) {                                                       
        // @todo : this is the world's worst way to increment                      
        // the number of equations we are solving. Let's fix                       
        // this before we embarrass ourselves                                      
        this->incrementNumEquations();                                             
        this->setSolveSumForceEquationInDirectionY( index , true );                
      } else if ( this->GetNodeType( index ) == Connect ) {                          
        // The node displacement and applied for cannot                            
        // be simultaneously defined. Here we issue a warning                      
        // message                                                                 
        throw MAP_ERROR_30;                                                        
      }
    } else {                                                                         
      // The node does not meet any initialization criteria                        
      throw MAP_ERROR_31;                                                          
    }
  } catch ( MAP_ERROR_CODE &code ) {                                               
    std::string str = "";                                                          
    std::ostringstream S;                                                          
    S << index+1;                                                                  
    str += boost::lexical_cast < std::string > ( code );                           
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";           

    Err.set_error_key( MAP_ERROR );                                                

    // Here is the exception message                                               
    Msg.RecordToErrorList( str );                                                  

    // create context to write error in summary.map                                
    str.erase(2,1);                                                                
    S.str("");S.clear();                                                           
    S << ">>>> " +str +"    :    MAP error in file ";                              
    S << __FILE__;                                                                 
    S << " on line ";                                                              
    S << __LINE__;                                                                 
    S << " occurred.";                                                             
    Msg.WriteErrorToOutputFile( S.str() );                                         
  }

  // ==========   Z Node VarType assignment  ===============     <-------------------------------------------------------------+
  // Assign data in the Node                                                                             
  try {                                                                                                  
    // If Z has a fixed value and it IS NOT attached to a vessel,                                        
    // set it as a parameter                                                                             
    if ( this->nodeVarTypeBool( index , &Node::Z ) == true && this->GetNodeType( index ) != Vessel ) {   
      // assign as a MAP_ParameterType_class                                                             
      this->setNodeReferenceToMAPType( index , P, &Node::Z );                                                
    } else if ( this->nodeVarTypeBool( index , &Node::Z ) == true && this->GetNodeType( index ) == Vessel ) {  
      // If Z has a fixed value and it is attached to a vessel, set it as an input 
      // assign as a MAP_InputType_class                                                                     
      this->setNodeReferenceToMAPType( index , I , &Node::Z );                                               
    } else { // if Z is variable (iterated), then it must be a constraint                                             
      // @todo : delete this in the next version of MAP after we verify this                                 
      // logic is correct                                                                                    
      assert( this->GetNodeType( index ) == Connect );                                                       

      // ensure this node is a 'Connect'                                                                     
      if( this->GetNodeType( index ) != Connect ) {                                                          
        // @todo : remove assertion in the next release of MAP                                               
        // once we verify things are working correctly                                                       

        throw MAP_ERROR_27;                                                                                  
      }

      // assign as a MAP_ConstraintStateType_class                                                           
      this->setNodeReferenceToMAPType( index , C , &Node::Z );                                         

      // Increment the number of equations by                                                          
      // 1 only if the following conditions are                                                        
      // met:                                                                                          
      //   -- Z.is_fixed  = false                                                                      
      //   -- FZ.is_fixed = true                                                                       
      //   -- Node.type   = Connect                                                                    
      if ( this->nodeVarTypeBool( index, &Node::Z )==false && this->nodeVarTypeBool( index, &Node::FZ )==true 
           && this->GetNodeType( index ) == Connect ) {                               
        // @todo : this is the world's worst way to increment                         
        // the number of equations we are solving. Let's fix                          
        // this before we embarrass ourselves                                         
        this->incrementNumEquations();                                                
        this->setSolveSumForceEquationInDirectionZ( index , true );                   
      } else {
        // The node displacement and applied for cannot                               
        // be simultaneously defined. Here we issue an error                          
        // message                                                                    
        throw MAP_ERROR_32;                                                           
      }
    }
  } catch ( MAP_ERROR_CODE &code ) {                                                  
    std::string str = "";                                                             
    std::ostringstream S;                                                             
    S << index+1;                                                                     
    str += boost::lexical_cast < std::string > ( code );                              
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";              

    Err.set_error_key( MAP_ERROR );                                                   

    // Here is the exception message                                                  
    Msg.RecordToErrorList( str );                                                     

    // create context to write error in summary.map                                   
    str.erase(2,1);                                                                   
    S.str("");S.clear();                                                              
    S << ">>>> " +str +"    :    MAP error in file ";                                 
    S << __FILE__;                                                                    
    S << " on line ";                                                                 
    S << __LINE__;                                                                    
    S << " occurred.";                                                                
    Msg.WriteErrorToOutputFile( S.str() );                                            
  }

  // ==========   M Node VarType assignment  ===============     <-------------------------------------------------------------+
  // Assign data in the Node           
  try {                                                                                                     
    if ( this->nodeVarTypeBool( index , &Node::M ) == true ) {                                              
      // assign as a MAP_ParameterType_class                                                                
      this->setNodeReferenceToMAPType( index , P, &Node::M );                                              
    } else { 
      // assign as a MAP_ConstraintStateType_class                                                         
      this->setNodeReferenceToMAPType( index , C , &Node::M );                                             

      // Increment the number of equations by 1 only if the following conditions are met:                  
      //   -- M.is_fixed  = false                                                                          
      //   -- B.is_fixed  = true                                                                           
      //   -- Z.is_fixe d = true                                                                           
      //   -- FZ.is_fixed = true                                                                           
      //   -- Node.type   = Connect                                                                        
      if ( this->nodeVarTypeBool( index, &Node::M )==false && this->nodeVarTypeBool( index, &Node::Z )==true 
           && this->nodeVarTypeBool( index, &Node::B )==true && this->nodeVarTypeBool( index, &Node::FZ )==true
           && this->GetNodeType( index ) == Connect ) {                                                        
        this->incrementNumEquations();                                                          
      }                                                                                         
      else if ( this->GetNodeType( index ) == Connect ) {                                       
        // The node displacement and applied for cannot be simultaneously defined. Here we      
        // issue a warning message                                                              
        throw MAP_ERROR_33;                                                                     
      }

      if ( this->GetNodeType( index ) != Fix ) {                                                
        throw MAP_ERROR_35;                                                                     
      }
    }
  } catch ( MAP_ERROR_CODE &code ) {                                                            
    std::string str = "";                                                                       
    std::ostringstream S;                                                                       
    S << index+1;                                                                               
    str += boost::lexical_cast < std::string > ( code );                                        
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";                        

    Err.set_error_key( MAP_ERROR );                                                             

    // Here is the exception message                                                            
    Msg.RecordToErrorList( str );                                                               

    // create context to write error in summary.map          
    str.erase(2,1);                                          
    S.str("");S.clear();                                     
    S << ">>>> " +str +"    :    MAP error in file ";        
    S << __FILE__;                                           
    S << " on line ";                                        
    S << __LINE__;                                           
    S << " occurred.";                                       
    Msg.WriteErrorToOutputFile( S.str() );                   
  }

  // ==========   B Node VarType assignment  ===============     <-------------------------------------------------------------+
  // Assign data in the Node   
  try {                                                                
    if ( this->nodeVarTypeBool( index , &Node::B ) == true ) {         
      // assign as a MAP_ParameterType_class                           
      this->setNodeReferenceToMAPType( index , P, &Node::B );          
    } else {                                                             
      // assign as a MAP_ConstraintStateType_class                     
      this->setNodeReferenceToMAPType( index , C , &Node::B );         

      // Increment the number of equations by 1 only if the following conditions are met:                      
      //   -- M.is_fixed  = true                                                                               
      //   -- B.is_fixed  = false                                                                              
      //   -- Z.is_fixed  = true                                                                               
      //   -- FZ.is_fixed = true                                                                               
      //   -- Node.type   = Connect                                                                            
      if ( this->nodeVarTypeBool( index , &Node::B )==false && this->nodeVarTypeBool( index , &Node::Z )==true 
           && this->nodeVarTypeBool( index, &Node::M )==true && this->nodeVarTypeBool( index, &Node::FZ )==true
           && this->GetNodeType( index ) == Connect ) {                                                        
        this->incrementNumEquations();                                                                       
      }                                                                                                      
      else if ( this->GetNodeType( index ) == Connect ) {                                                    
        // The node displacement and applied for cannot be simultaneously defined. Here we                   
        // issue a warning message                                                                           
        throw MAP_ERROR_34;                                                                                  
      }

      if ( this->GetNodeType( index ) != Fix ) {                                                             
        throw MAP_ERROR_35;                                                                                  
      }
    }
  } catch ( MAP_ERROR_CODE &code ) {                                                                         
    std::string str = "";                                                                                    
    std::ostringstream S;                                                                                    
    S << index+1;                                                                                            
    str += boost::lexical_cast < std::string > ( code );                                                     
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";                                     

    Err.set_error_key( MAP_ERROR );                                                                          

    // Here is the exception message                                                                         
    Msg.RecordToErrorList( str );                                                                            

    // create context to write error in summary.map                                                          
    str.erase(2,1);                                                                                          
    S.str("");S.clear();                                                                                     
    S << ">>>> " +str +"    :    MAP error in file ";                                                        
    S << __FILE__;                                                                                           
    S << " on line ";                                                                                        
    S << __LINE__;                                                                                           
    S << " occurred.";                                                                                       
    Msg.WriteErrorToOutputFile( S.str() );                                                                   
  }

  /**
   * ===============  FX, FY and FZ rules  ======================================
   * 	   
   * This part is a little tricky. FX, FY and FZ can be set as iterated values in 
   * the MAP input file; however, they are not actually iterated (solved). If FX, 
   * FY or FZ is selected as an iterated value, it is a relfection of the sum forces 
   * applied to the node, H and V (once converted into the inertial/global frame).
   * 
   * The process below sets H and V at the element level as iterated variables.
   * =============================================================================
   */

  // ==========   FX Node VarType assignment  ===============     <------------------------------------------------------------+
  // Assign data in the Node. FX can only be a parameter or constraint                                     
  try { 
    if ( this->nodeVarTypeBool( index , &Node::FX ) == true ) {                                            
      // assign as a MAP_ParameterType_class                                                               
      this->setNodeReferenceToMAPType( index , P, &Node::FX );                                             
    } else { 
      setNodeReferenceToMAPType( index , OtherStateDataTypes , &Node::FX );

      // Increment the number of equations by 1 only if the following conditions are met:                      
      //   -- FX.is_fixed = false                                                                              
      //   -- X.is_fixed  = true                                                                               
      //   -- Node.type   = Connect                                                                            
      if ( nodeVarTypeBool( index, &Node::FX )==false && nodeVarTypeBool( index, &Node::X )==true  
           && GetNodeType( index ) == Connect ) {                                                        
        this->incrementNumEquations();                                                                       
      } else if ( this->GetNodeType( index ) == Connect ) {                                                      
        throw MAP_ERROR_36;                                                                                  
      }
    }

    // if the node is attached to a vessel, then we must output                                                
    // the node force                                                                                          
    if ( this->GetNodeType( index ) == Vessel ){                                                               
      this->setNodeReferenceToMAPType( index , O , &Node::FX );                                                
    }
  } catch ( MAP_ERROR_CODE &code ) {                                      
    std::string str = "";                                                 
    std::ostringstream S;                                                 
    S << index+1;                                                         
    str += boost::lexical_cast < std::string > ( code );                  
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";  
 
    Err.set_error_key( MAP_ERROR );                                       
 
    // Here is the exception message                                      
    Msg.RecordToErrorList( str );                                         
 
    // create context to write error in summary.map                       
    str.erase(2,1);                                                       
    S.str("");S.clear();                                                  
    S << ">>>> " +str +"    :    MAP error in file ";                     
    S << __FILE__;                                                        
    S << " on line ";                                                     
    S << __LINE__;                                                        
    S << " occurred.";                                                    
    Msg.WriteErrorToOutputFile( S.str() );                                
  }

  // ==========   FY Node VarType assignment  ===============     <------------------------------------------------------------+
  // Assign data in the Node. FY can only be a parameter or constraint 
  try {                                                                                  
    if ( this->nodeVarTypeBool( index , &Node::FY ) == true ) {                          
      // assign as a MAP_ParameterType_class                                             
      this->setNodeReferenceToMAPType( index , P, &Node::FY );                           
    } else {                                                                               
      this->setNodeReferenceToMAPType( index , this->OtherStateDataTypes , &Node::FY );

      // Increment the number of equations by 1 only if the following conditions are met:                    
      //   -- FY.is_fixed = false                                                                            
      //   -- Y.is_fixed  = true                                                                             
      //   -- Node.type   = Connect                                                                          
      if ( this->nodeVarTypeBool( index, &Node::FY )==false && this->nodeVarTypeBool( index, &Node::Y )==true
           && this->GetNodeType( index ) == Connect ) {                                                      
        this->incrementNumEquations();                                                                       
      }                                                                                                      
      else if ( this->GetNodeType( index ) == Connect ) {                                                    
        throw MAP_ERROR_37;                                                                                  
      }
    }

    // if the node is attached to a vessel, then we must output                                              
    // the node force                                                                                        
    if ( this->GetNodeType( index ) == Vessel ){                                                             
      this->setNodeReferenceToMAPType( index , O , &Node::FY );                                              
    }
  } catch ( MAP_ERROR_CODE &code ) {                                                                         
    std::string str = "";                                                                                    
    std::ostringstream S;                                                                                    
    S << index+1;                                                                                            
    str += boost::lexical_cast < std::string > ( code );                                                     
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";                                     

    Err.set_error_key( MAP_ERROR );                                                                          

    // Here is the exception message                                                                         
    Msg.RecordToErrorList( str );                                                                            

    // create context to write error in summary.map                                                          
    str.erase(2,1);                                                                                          
    S.str("");S.clear();                                                                                     
    S << ">>>> " +str +"    :    MAP error in file ";                                                        
    S << __FILE__;                                                                                           
    S << " on line ";                                                                                        
    S << __LINE__;                                                                                           
    S << " occurred.";                                                                                       
    Msg.WriteErrorToOutputFile( S.str() );                                                                   
  }

  // ==========   FZ Node VarType assignment  ===============     <------------------------------------------------------------+
  // Assign data in the Node. FZ can only be a parameter or constraint
  try {                                                                                             
    if ( this->nodeVarTypeBool( index , &Node::FZ ) == true ) {                                     
      // assign as a MAP_ParameterType_class                                                        
      this->setNodeReferenceToMAPType( index , P, &Node::FZ );                                      
    } else {                                                                                          
      this->setNodeReferenceToMAPType( index , this->OtherStateDataTypes , &Node::FZ );

      // Increment the number of equations by 1 only if the following conditions are met: 
      //   -- FZ.is_fixed = false                                                         
      //   -- Z.is_fixed  = true                                                          
      //   -- B.is_fixed  = true                                                          
      //   -- M.is_fixed  = true                                                          
      //   -- Node.type   = Connect                                                       
      if ( this->nodeVarTypeBool( index, &Node::FZ )==false && this->nodeVarTypeBool( index, &Node::B )==true  
           && this->nodeVarTypeBool( index, &Node::M )==true && this->nodeVarTypeBool( index, &Node::Z )==true 
           && this->GetNodeType( index ) == Connect ) {                                                        
        this->incrementNumEquations();                                                                       
      } else if ( this->GetNodeType( index ) == Connect ) {                                                      
        throw MAP_ERROR_38;                                                                                  
      };                                                                                                       
    };                                                                                                         

    // if the node is attached to a vessel, then we must output                                                
    // the node force                                                                                          
    if ( this->GetNodeType( index ) == Vessel ){                                                               
      this->setNodeReferenceToMAPType( index , O , &Node::FZ );                                                
    }
  } catch ( MAP_ERROR_CODE &code ) {                                                                           
    std::string str = "";                                                                                      
    std::ostringstream S;                                                                                      
    S << index+1;                                                                                              
    str += boost::lexical_cast < std::string > ( code );                                                       
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";                                       

    Err.set_error_key( MAP_ERROR );                                                                            

    // Here is the exception message                                                                           
    Msg.RecordToErrorList( str );                                                                              

    // create context to write error in summary.map                             
    str.erase(2,1);                                                             
    S.str("");S.clear();                                                        
    S << ">>>> " +str +"    :    MAP error in file ";                           
    S << __FILE__;                                                              
    S << " on line ";                                                           
    S << __LINE__;                                                              
    S << " occurred.";                                                          
    Msg.WriteErrorToOutputFile( S.str() );                                      
  }

  // final check to make sure we are not defining too many variables being solved (more equations than unknowns)
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
    }
  } catch( MAP_ERROR_CODE &code ) {
    std::string str = "";
    std::ostringstream S;
    S << index+1;
    str += boost::lexical_cast < std::string > ( code );
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";

    Err.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordToErrorList( str );
           
    // create context to write error in summary.map 
    str.erase(2,1);                                
    S.str("");S.clear();                            
    S << ">>>> " +str +"    :    MAP error in file ";
    S << __FILE__;                                  
    S << " on line ";                               
    S << __LINE__;                                  
    S << " occurred.";                              
    Msg.WriteErrorToOutputFile( S.str() );                     
  }

  // Make sure X direction doesn't have more unknown than equations
  counter = 0;
  try {
    if ( this->nodeVarTypeBool( index , &Node::X  ) == false ) counter++;
    if ( this->nodeVarTypeBool( index , &Node::FX ) == false ) counter++;

    if (counter > 1) {
      // The node displacement and applied for cannot  be simultaneously defined. Here we issue a warning message
      throw MAP_ERROR_40;
    }
  } catch( MAP_ERROR_CODE &code ) {
    std::string str = "";
    std::ostringstream S;
    S << index+1;
    str += boost::lexical_cast < std::string > ( code );
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";

    Err.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordToErrorList( str );

    // create context to write error in summary.map 
    str.erase(2,1);                                
    S.str("");S.clear();                            
    S << ">>>> " +str +"    :    MAP error in file ";
    S << __FILE__;                                  
    S << " on line ";                               
    S << __LINE__;                                  
    S << " occurred.";                              
    Msg.WriteErrorToOutputFile( S.str() );                     
  }

  // Make sure Y direction doesn't have more unknown than equations
  counter = 0;
  try {
    if ( this->nodeVarTypeBool( index , &Node::Y  ) == false ) counter++;
    if ( this->nodeVarTypeBool( index , &Node::FY ) == false ) counter++;

    if (counter > 1){
      // The node displacement and applied for cannot be simultaneously defined. Here we issue a warning  
      // message            
      throw MAP_ERROR_41;
    }
  } catch( MAP_ERROR_CODE &code ) {
    std::string str = "";
    std::ostringstream S;
    S << index+1;
    str += boost::lexical_cast < std::string > ( code );
    str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";

    Err.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordToErrorList( str );

    // create context to write error in summary.map 
    str.erase(2,1);                                
    S.str("");S.clear();                            
    S << ">>>> " +str +"    :    MAP error in file ";
    S << __FILE__;                                  
    S << " on line ";                               
    S << __LINE__;                                  
    S << " occurred.";                              
    Msg.WriteErrorToOutputFile( S.str() );                     
  }
};



// ============================================================================
// GetSolverOptionsString( )
// ============================================================================
std::string MAP_OtherStateType_class::
GetSolverOptionsString( MAP_InitInputType_class &Init )
{
  std::ostringstream S;

  S << "MAP Solver options:\n";

  for ( int i=0 ; i<Init.SizeofSolverOptions() ; i++ ){
    S << "    " + Init.GetSolverOptions(i);
  }
  return S.str();
}


// ============================================================================
// SetVar
//
// only used to set 'Connect' node FX, FY, and FZ 
// ============================================================================
void MAP_OtherStateType_class::
SetVar( const int index , 
        double    value )
{
  OtherStateDataTypes.SetVar( index , value );
} 


// ============================================================================
// writeLinearizedStiffnessMatrix
//
// Computes the linearized stiffness matrix. This result will be writtent to
// the MAP output file.
//
// The linearization model used centered finite differencing:
//  
//   -- [K]  =  [ F(X+epsilon) - F(X-epsilon) ]/(2*epsilon)
//
// @todo : In the future, the value for epsilon will be a run-time command. For
//         now, it is set at compile time. 
// @todo : In the cross-product term, check if the the element Hx, Hy and V 
//         force should be used, or if the node force FX, FY and FZ should be 
//         used. At the moment, FX, FY and HZ are being used. 
// ============================================================================
void MAP_OtherStateType_class::
writeLinearizedStiffnessMatrix( MAP_Message_class &Msg )
{

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

  MAP_ErrStat_class        Err; 
  std::ostringstream S;

  S << std::fixed << std::setprecision(4);

  for ( unsigned int i=0 ; i<node.size() ; i++ ) {                              
    if( node[i]->type == Vessel ){			      
      // For Vessel fairleads: if FX or FY value does not begin with a '#', 
      if ( node[i]->GetVarTypeFixedBool( &Node::X )  == false ||
           node[i]->GetVarTypeFixedBool( &Node::Y )  == false ||
           node[i]->GetVarTypeFixedBool( &Node::Z )  == false ||
           node[i]->GetVarTypeFixedBool( &Node::FX ) == true  ||
           node[i]->GetVarTypeFixedBool( &Node::FY ) == true  ||
           node[i]->GetVarTypeFixedBool( &Node::FZ ) == true   )
      {
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
  if (flg) {
    for ( int jj=0 ; jj<SIX ; jj++ ) {

      /**
       * =======  Foward finite difference:  ======     <----------------------------------------------------------+
       *                                                                                    
       *  F( X + epsilon )                                                                  
       */                                                                                   

      for ( unsigned int i=0 ; i<node.size() ; i++ ) {                                      
        if ( node[i]->type == Vessel ) {						    

          // Store the original vessel displacement in some varaible. This will be restored 
          // once the new displacement+epsilon is solved. 				    

          // perturb the vessel by epsilon						    
          if ( jj == 0 ) { // X + epsilon						    
            delta = node[i]->GetPriorXValue( )+epsilon;					 
            node[i]->SetVarTypeValue( &Node::X , delta );				    
          }										    
          else if ( jj == 1 ) { // Y + epsilon						    
            delta = node[i]->GetPriorYValue( )+epsilon;					 
            node[i]->SetVarTypeValue( &Node::Y , delta );				    
          } else if ( jj == 2 ) { // Z + epsilon						    
            delta = node[i]->GetPriorZValue( )+epsilon;					 
            node[i]->SetVarTypeValue( &Node::Z , delta );				    
          } else if ( jj == 3 ) { // phi + epsilon					    

            // -(phi*z) about j axis 							    
            //  (phi*y) about k axis							    
            delta   = node[i]->GetPriorYValue( ) - node[i]->GetPriorZValue( )*epsilon; 					 
            delta_2 = node[i]->GetPriorZValue( ) + node[i]->GetPriorYValue( )*epsilon; 					

            node[i]->SetVarTypeValue( &Node::Y , delta   ); 				    
            node[i]->SetVarTypeValue( &Node::Z , delta_2 ); 				    
          } else if ( jj == 4 ) { // theta + epsilon					    

            //  (theta*z) about i axis 							    
            // -(theta*x) about k axis							    
            delta   = node[i]->GetPriorXValue( ) + node[i]->GetPriorZValue( )*epsilon; 	
            delta_2 = node[i]->GetPriorZValue( ) - node[i]->GetPriorXValue( )*epsilon; 	

            node[i]->SetVarTypeValue( &Node::X , delta   ); 				    
            node[i]->SetVarTypeValue( &Node::Z , delta_2 ); 				    
          } else if ( jj == 5 ) { // psi + epsilon					    

            // -(psi*y) about i axis 							    
            //  (psi*x) about j axis							    
            delta   = node[i]->GetPriorXValue( ) - node[i]->GetPriorYValue( )*epsilon; 	
            delta_2 = node[i]->GetPriorYValue( ) + node[i]->GetPriorXValue( )*epsilon; 	

            node[i]->SetVarTypeValue( &Node::X , delta   ); 
            node[i]->SetVarTypeValue( &Node::Y , delta_2 ); 
          }
        }
      }
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
      }

      /**
       * =======  Mooring force with Forward finite difference  ======     <---------------------------------------+
       *                                                                
       * cross product to calculate moments :                           
       *								
       *   | i   j   k  |      Y*Fz  - Z*Fy				
       *   | X   Y   Z  |  =   Z*Fx  - X*Fz				
       *   | Fx  Fy  Fz |      X*Fy  - Y*Fx				
       */								
									
      for ( unsigned int i=0 ; i<node.size() ; i++ ) {			
        if ( node[i]->type == Vessel ) {				
          force[0] += node[i]->GetVarTypeValue( &Node::FX );  		
          force[1] += node[i]->GetVarTypeValue( &Node::FY );		
          force[2] += node[i]->GetVarTypeValue( &Node::FZ );		
      
          Mx = ((node[i]->GetVarTypeValue( &Node::FZ ) ) *(node[i]->GetVarTypeValue( &Node::Y ))	
                - (node[i]->GetVarTypeValue( &Node::FY ))*(node[i]->GetVarTypeValue( &Node::Z )));
          My = ((node[i]->GetVarTypeValue( &Node::FX ))  *(node[i]->GetVarTypeValue( &Node::Z ))	
                - (node[i]->GetVarTypeValue( &Node::FZ ))*(node[i]->GetVarTypeValue( &Node::X )));
          Mz = ((node[i]->GetVarTypeValue( &Node::FY ))  *(node[i]->GetVarTypeValue( &Node::X ))	
                - (node[i]->GetVarTypeValue( &Node::FX ))*(node[i]->GetVarTypeValue( &Node::Y )));	

          if ( jj == 0 || jj == 1 || jj == 2 ) {
            Mx = ((node[i]->GetVarTypeValue( &Node::FZ ))  *(node[i]->GetPriorYValue( ))	   
                  - (node[i]->GetVarTypeValue( &Node::FY ))*(node[i]->GetPriorZValue( )));  
            My = ((node[i]->GetVarTypeValue( &Node::FX ))  *(node[i]->GetPriorZValue( ))	   
                  - (node[i]->GetVarTypeValue( &Node::FZ ))*(node[i]->GetPriorXValue( )));  
            Mz = ((node[i]->GetVarTypeValue( &Node::FY ))  *(node[i]->GetPriorXValue( ))	   
                  - (node[i]->GetVarTypeValue( &Node::FX ))*(node[i]->GetPriorYValue( )));  

            force[3] += Mx;
            force[4] += My;
            force[5] += Mz;
          } else if ( jj == 3 ) {
            force[3] +=  Mx;
            force[4] +=  My*cos(epsilon) - Mz*sin(epsilon);	
            force[5] +=  My*sin(epsilon) + Mz*cos(epsilon);
          } else if ( jj == 4) {
            force[3] +=  Mx*cos(epsilon) + Mz*sin(epsilon);
            force[4] +=  My;	
            force[5] +=  -Mx*sin(epsilon) + Mz*cos(epsilon);
          } else if ( jj == 5 ) {
            force[3] +=  Mx*sin(epsilon) + My*cos(epsilon);	
            force[4] +=  Mx*cos(epsilon) - My*sin(epsilon);	
            force[5] +=  Mz;
          }

          Mx = 0.0;
          My = 0.0;
          Mz = 0.0;
        }
      }
      //============== <END> =======================================================================================
								
      
      // Restore the original position of the vessel nodes prior to perturbing it by epsilon					
      for ( unsigned int i=0 ; i<node.size() ; i++ ){				
        if ( node[i]->type == Vessel ){					
          node[i]->SetVarTypeValue( &Node::X , node[i]->GetPriorXValue() );	
          node[i]->SetVarTypeValue( &Node::Y , node[i]->GetPriorYValue() );	
          node[i]->SetVarTypeValue( &Node::Z , node[i]->GetPriorZValue() );	
        }
      }

      // =======  Backward finite difference:  ======     <--------------------------------------------------------+
      //  F( X - epsilon )
      for ( unsigned int i=0 ; i<node.size() ; i++ ) {  
        // Only Vessel nodes are moded, because we want the stiffness to be relative to small vessel displacements
        if ( node[i]->type == Vessel ) {
          // Store the original vessel displacement in some varaible. This will be restored
          // once the new displacement+epsilon is solved.

          // perturb the vessel by epsilon
          if ( jj == 0 ) { // X - epsilon
            delta = node[i]->GetPriorXValue( )-epsilon;	
            node[i]->SetVarTypeValue( &Node::X , delta );		
          } else if ( jj == 1 ) { // Y - epsilon
            delta = node[i]->GetPriorYValue( )-epsilon;
            node[i]->SetVarTypeValue( &Node::Y , delta );		
          } else if ( jj == 2 ) { // Z - epsilon			
            delta = node[i]->GetPriorZValue( )-epsilon;	
            node[i]->SetVarTypeValue( &Node::Z , delta );		
          } else if ( jj == 3 ) { // phi  - epsilon			

            //  (phi*z) about j axis 					
            // -(phi*y) about k axis					
            delta   = node[i]->GetPriorYValue( ) 			
              + node[i]->GetPriorZValue( )*epsilon; 			
            delta_2 = node[i]->GetPriorZValue( ) 			
              - node[i]->GetPriorYValue( )*epsilon; 		

            node[i]->SetVarTypeValue( &Node::Y , delta   ); 	
            node[i]->SetVarTypeValue( &Node::Z , delta_2 ); 	
          } else if ( jj == 4 ) { // theta - epsilon		

            // -(theta*z) about i axis 				
            //  (theta*x) about k axis				
            delta   = node[i]->GetPriorXValue( ) 		
              - node[i]->GetPriorZValue( )*epsilon; 		
            delta_2 = node[i]->GetPriorZValue( ) 		
              + node[i]->GetPriorXValue( )*epsilon; 		

            node[i]->SetVarTypeValue( &Node::X , delta   ); 	
            node[i]->SetVarTypeValue( &Node::Z , delta_2 ); 	
          } else if ( jj == 5 ) { // psi - epsilon		

            //  (psi*y) about i axis 				
            // -(psi*x) about j axis				
            delta   = node[i]->GetPriorXValue( ) 		
              + node[i]->GetPriorYValue( )*epsilon; 		
            delta_2 = node[i]->GetPriorYValue( ) 		
              - node[i]->GetPriorXValue( )*epsilon; 		
            node[i]->SetVarTypeValue( &Node::X , delta   ); 	
            node[i]->SetVarTypeValue( &Node::Y , delta_2 ); 	
          }
        }
      }
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
      }
	
	

      //  =======  Mooring force with Backward finite difference  ======     <--------------------------------------+
      //                                                         
      //  cross product to calculate moments :                   
      // 							
      //    | i   j   k  |      Y*Fz  - Z*Fy			
      //    | X   Y   Z  |  =   Z*Fx  - X*Fz			
      //    | Fx  Fy  Fz |      X*Fy  - Y*Fx			
      for ( unsigned int i=0 ; i<node.size() ; i++ ) {	
        if ( node[i]->type == Vessel ) {			
          force[0] -= node[i]->GetVarTypeValue( &Node::FX );	
          force[1] -= node[i]->GetVarTypeValue( &Node::FY );	
          force[2] -= node[i]->GetVarTypeValue( &Node::FZ );	

          Mx = ((node[i]->GetVarTypeValue( &Node::FZ ) ) *(node[i]->GetVarTypeValue( &Node::Y ))	
                - (node[i]->GetVarTypeValue( &Node::FY ))*(node[i]->GetVarTypeValue( &Node::Z )));
          My = ((node[i]->GetVarTypeValue( &Node::FX ))  *(node[i]->GetVarTypeValue( &Node::Z ))	
                - (node[i]->GetVarTypeValue( &Node::FZ ))*(node[i]->GetVarTypeValue( &Node::X )));
          Mz = ((node[i]->GetVarTypeValue( &Node::FY ))  *(node[i]->GetVarTypeValue( &Node::X ))	
                - (node[i]->GetVarTypeValue( &Node::FX ))*(node[i]->GetVarTypeValue( &Node::Y )));	

          if ( jj == 0 || jj == 1 || jj == 2 ){
            Mx = ((node[i]->GetVarTypeValue( &Node::FZ ) ) *(node[i]->GetPriorYValue( ))	
                  - (node[i]->GetVarTypeValue( &Node::FY ))*(node[i]->GetPriorZValue( )));	
            My = ((node[i]->GetVarTypeValue( &Node::FX ))  *(node[i]->GetPriorZValue( ))	
                  - (node[i]->GetVarTypeValue( &Node::FZ ))*(node[i]->GetPriorXValue( )));	
            Mz = ((node[i]->GetVarTypeValue( &Node::FY ))  *(node[i]->GetPriorXValue( ))
                  - (node[i]->GetVarTypeValue( &Node::FX ))*(node[i]->GetPriorYValue( )));	

            force[3] -= Mx;
            force[4] -= My;
            force[5] -= Mz;
          } else if ( jj == 3 ) {
            force[3] -=  Mx;
            force[4] -=  My*cos(epsilon) - Mz*sin(epsilon);	
            force[5] -=  My*sin(epsilon) + Mz*cos(epsilon);
          } else if ( jj == 4) {
            force[3] -=  Mx*cos(epsilon) + Mz*sin(epsilon);
            force[4] -=  My;	
            force[5] -= -Mx*sin(epsilon) + Mz*cos(epsilon);
          } else if ( jj == 5 ) {
            force[3] -=  Mx*sin(epsilon) + My*cos(epsilon);	
            force[4] -=  Mx*cos(epsilon) - My*sin(epsilon);	
            force[5] -=  Mz;
          }

          Mx = 0.0;
          My = 0.0;
          Mz = 0.0;
        }
      }
      //============== <END> =======================================================================================
      
      // Restore the original position of the vessel nodes prior to perturbing it by epsilon
      for ( unsigned int i=0 ; i<node.size() ; i++ ) {
        if ( node[i]->type == Vessel ) {
          node[i]->SetVarTypeValue( &Node::X , node[i]->GetPriorXValue() );
          node[i]->SetVarTypeValue( &Node::Y , node[i]->GetPriorYValue() );
          node[i]->SetVarTypeValue( &Node::Z , node[i]->GetPriorZValue() );
        }
      }


      // restore original solution before epsilon displacement.
      //
      // 1) call numeric_method->PetscSolve( Error , Msg );
      // 
      // 2) reset the sum forces FX, FY and FZ in each node
      Solve( Err, Msg );  
      for( int i=0 ; i<getSizeOfNode() ; i++) {
        initializeCableNode( i );
      }

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

      assert( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      do {							      //             |
        write += " ";					      //             |
      } while( write.size()<_TEXT_COLOR::STR_LEN );		      //             |
      write_line_5 += write;                                        //             |
      write.clear();                                                //   ----------+
      //============================================================================

	
      /**
       * ==========   Yaw row of matrix [K]   ================     <---------------+
       * Prints the Yaw row linearized stiffness 		   
       *							   
       * To be sure the length of 'write_line_1' is not greater    
       * than length we are allow to across, this assert statement 
       * is raised.						   
       *							   
       * Fill spaces in the 'write' string with white spaces to    
       * occupy area in the output string (to make the output	   
       * print in alignment)					   
       */							   
								   
      S << force[5]/(2*epsilon);                                   
      if ( force[5]/(2*epsilon) >= 0.0) write += " ";		   
      write += S.str();
      S.str("");S.clear(); 					   

      assert( write.size()<_TEXT_COLOR::STR_LEN );		   
      do {							   
        write += " ";
      } while( write.size()<_TEXT_COLOR::STR_LEN );		   
      write_line_6 += write;
      write.clear(); 
      //============================================================================

      for( int i=0 ; i<SIX ; i++) {
        force[i] = 0.0;
      }
    }

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
  } else {
    Msg.WriteDataToOutputFile( "Linearized Stiffness Matrix\n" );
    Msg.WriteDataToOutputFile( "--------------------------------------------------------------------------" );
    Msg.WriteDataToOutputFile( "--------------------------------------------------------------------------\n" );
    Msg.WriteDataToOutputFile( "    Stiffness matrix could not be computed with given\n" );
    Msg.WriteDataToOutputFile( "    input file options\n" );
  }
};
