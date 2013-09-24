/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        MAP_InitInputType.cpp
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


#include "MAP_InitInputType.h"  // Preprocessor Defitions in Prerequisite.h:
                                //   #include <string>
                                //   #include <map>


// ============================================================================
// Set_Cable_Library_Data
//
// ============================================================================
void 
MAP_InitInputType_class::
SetCableLibraryData( const std::string &text ) 
{ 
  cable_library_data.push_back( text ); 
};


// ============================================================================
// Set_Solver_Options
//
// ============================================================================
void MAP_InitInputType_class::
SetSolverOptions( const std::string &text ) 
{ 
  solver_options.push_back( text ); 
};


// ============================================================================
// SetNodeData
//
// ============================================================================
void MAP_InitInputType_class::
SetNodeData( const std::string &text ) 
{ 
  node_data.push_back( text );
  this->size_of_nodes = node_data.size();
};


// ============================================================================
// SetCoupledToFastFlag
//
// ============================================================================
void MAP_InitInputType_class::
SetCoupledToFastFlag( const bool flag )
{
  this->FAST_coupled = flag;
};


// ============================================================================
// GetCoupledToFastFlag
//
// ============================================================================
bool MAP_InitInputType_class::
GetCoupledToFastFlag( ) const 
{
  return this->FAST_coupled;
};


// ============================================================================
// SetElementData 
//
// Uses the input line in the 'LINE PROPERTIERS' section of the MAP input file and
//  parses the line into individual words (using boost::split). 
//  
// @input : T  -- a single line in the 'LINE PROPERTIES' section of the MAP
//                input file.
//
// @todo  : Fix the exception in the else statement to handle cases where 'Repeat'
//          is mispelled, so a segfault does not occur. 
// ============================================================================
void MAP_InitInputType_class::
SetElementData( const std::string &text) 
{ 
  std::vector<std::string> repeat_line;
  
  //repeat_line[0]   -- element number
  //repeat_line[1]   -- line type (such as steel, nylon, ect.)
  //repeat_line[2]   -- unstretched element length
  //repeat_line[3]   -- anchor node
  //repeat_line[4]   -- fairlead node
  //repeat_line[5-N] -- option flags
  boost::split( repeat_line              , 
                text                     , 
                boost::is_any_of(" ")    , 
                boost::token_compress_on );
  
  // If the MAP input file line starts with the word 'repeat', then copy all the element 
  // previous to it. Otherwise, push the string back into the element so the properties
  // are stored
  if ( boost::iequals( repeat_line[1] , "Repeat")== false ) {
    element_data.push_back( text );      
    this->size_of_elements = element_data.size();
  } else if ( boost::iequals( repeat_line[1] , "Repeat")== true ) {

    // trim string and remove '\n' character at end of degree value. 
    repeat_line[2].erase( repeat_line[2].find_last_not_of(" \n")+1 );    

    double angle = boost::lexical_cast<double>( repeat_line[2] )*( MAP_CONST::PI/180 );

    std::string out = "";
    std::vector<std::string> words;
    std::ostringstream S;
    S << std::fixed << std::setprecision(3);
	
    std::string temp_X  = "";
    std::string temp_Y  = "";
    std::string temp_Z  = "";
    std::string temp_FX = "";
    std::string temp_FY = "";
    std::string temp_FZ = "";

    std::string X_string  = "";
    std::string Y_string  = "";
    std::string Z_string  = "";
    std::string FX_string = "";
    std::string FY_string = "";
    std::string FZ_string = "";

    bool X_is_fixed  = true;
    bool Y_is_fixed  = true;
    //bool Z_is_fixed  = true;  // @todo: figure out why this isn't needed
    bool FX_is_fixed = true;
    bool FY_is_fixed = true;
    //bool FZ_is_fixed = true;  // @todo: figure out why this isn't needed
	
    double X_double  = 0;
    double Y_double  = 0;
    double FX_double = 0;
    double FY_double = 0;

    int upper = 0;
    int lower = 0;

    for( int i=0 ; i<this->size_of_nodes ; i++){
      words.clear();

      // split words using boost library 
      boost::split( words                    , 
                    this->GetNodeData( i )   , 
                    boost::is_any_of(" \n")  , 
                    boost::token_compress_on );
	    
      S << size_of_nodes + node_counter;

      // node number
      out += S.str();
      out += " ";
      S.str("");S.clear();

      // node type
      out += words[1];
      out += " ";
      S.str("");S.clear();

      // X position
      if ( words[2][0] == '#') {
        X_is_fixed= false;
        // initial guess flag for MAP
        if ( boost::iequals( words[2] , "#") ) {
          temp_X = "#";
        } else {
          unsigned int i = 1;

          // Now cycle through the string and write it to the This parameter,
          // but ignore the firest character (which will be a '#')
          while(i < words[2].size() ){
            temp_X += words[2][i];
            i++;
          }
        }
      } else {
        temp_X = words[2];
      }

      // Y position
      if ( words[3][0] == '#') {
        Y_is_fixed= false;
        // initial guess flag for MAP
        if ( boost::iequals( words[3] , "#") ) temp_Y = "#";
        else {
          unsigned int i = 1;

          // Now cycle through the string and write it to the This parameter,
          // but ignore the firest character (which will be a '#')
          while(i < words[3].size() ){
            temp_Y += words[3][i];
            i++;
          }
        }
      } else {
        temp_Y = words[3];
      }

      X_double = boost::lexical_cast<double>( temp_X );
      Y_double = boost::lexical_cast<double>( temp_Y );

      // this is used in the event "repeat angle" is called in the MAP input file
      temp_X = boost::lexical_cast<std::string>(  X_double*cos(angle) + Y_double*sin(angle) );
      temp_Y = boost::lexical_cast<std::string>( -X_double*sin(angle) + Y_double*cos(angle) );

      // append '#' symbol to numeric value if variable is iterated
      if( X_is_fixed == false) {
        S << "#" + temp_X;
      } else {
        S << temp_X;
      }

      out += S.str();
      out += " ";	   
      S.str("");S.clear(); 

      // append '#' symbol to numeric value if variable is iterated
      if( Y_is_fixed == false) S << "#" + temp_Y;
      else S << temp_Y;
      out += S.str();
      out += " ";	   
      S.str("");S.clear(); 

      // Z position
      out += words[4];
      out += " ";	   
      S.str("");S.clear(); 

      // M point mass
      out += words[5];
      out += " ";	   
      S.str("");S.clear(); 

      // B buoyancy module
      out += words[6];
      out += " ";	   
      S.str("");S.clear(); 

      // FX force
      if ( words[7][0] == '#') {
        FX_is_fixed= false;

        // initial guess flag for MAP
        if ( boost::iequals( words[7] , "#") ) temp_FX = "#";
        else {
          unsigned int i = 1;

          // Now cycle through the string and write it to the This parameter,
          // but ignore the firest character (which will be a '#')
          while(i < words[7].size() ){
            temp_FX += words[7][i];
            i++;
          }
        }
      } else { 
        // FX is a fixed value (and does not start with '#')
        temp_FX = words[7];
      }
      
      // FY force
      if ( words[8][0] == '#' ) {
        FY_is_fixed= false;
        if ( boost::iequals( words[8] , "#") ){
          temp_FY = "#";
        } else {
          unsigned int i = 1;
          
          // Now cycle through the string and write it to the This parameter,
          //  but ignore the firest character (which will be a '#')
          while(i < words[8].size() ) {
            temp_FY += words[8][i];
            i++;
          }
        }
      } else {
        // FY is a fixed value (and does not start with '#')
        temp_FY = words[8];
      }

      try {
        FX_double = boost::lexical_cast<double>( temp_FX );
      } catch (boost::bad_lexical_cast const& ) {
        FX_double = 0;
      }

      try {
        FY_double = boost::lexical_cast<double>( temp_FY );
      } catch (boost::bad_lexical_cast const& ) {
        FY_double = 0;
      }

      if ( boost::iequals( words[7] , "#")==true && boost::iequals( words[8] , "#")==true ){
        temp_FX = "";
        temp_FY = "";
      } else {
        temp_FX = boost::lexical_cast<std::string>(  FX_double*cos(angle) + FY_double*sin(angle) );
        temp_FY = boost::lexical_cast<std::string>( -FX_double*sin(angle) + FY_double*cos(angle) );
        assert( boost::iequals( words[7] , "#")!=true && boost::iequals( words[8] , "#")!=true );
      }

      if(FX_is_fixed == false) {
        S << "#"+temp_FX;
      } else {
        S << temp_FX;
      }
      out += S.str();
      out += " ";	   
      S.str("");S.clear(); 

      if(FY_is_fixed == false) {
        S << "#"+temp_FY;
      } else {
        S << temp_FY;
      }
      out += S.str();
      out += " ";	   
      S.str("");S.clear(); 

      // FZ force
      out += words[9];
      out += " ";	   
      S.str("");S.clear(); 
	    
      temp_X.clear();
      temp_Y.clear();
      temp_Z.clear();
      temp_FX.clear();
      temp_FY.clear();
      temp_FZ.clear();
	    
      X_string.clear();
      Y_string.clear();
      Z_string.clear();
      FX_string.clear();
      FY_string.clear();
      FZ_string.clear();
	    
      X_is_fixed  = true;
      Y_is_fixed  = true;
      //Z_is_fixed  = true;
      FX_is_fixed = true;
      FY_is_fixed = true;
      //FZ_is_fixed = true;

      out +="\n";
      node_data.push_back( out );      
      out.clear();

      node_counter++;	    
    }

    // make sure [this->SizeofElementData()/this->size_of_elements] divides into a whole number
    assert( std::fmod( (this->SizeofElementData()/this->size_of_elements),1.0) == 0 );

    int I = this->SizeofElementData()/this->size_of_elements;

    for( int i=0 ; i<this->size_of_elements ; i++){
      words.clear();

      // split words using boost library
      boost::split( words                     , 
                    this->GetElementData( i ) , 
                    boost::is_any_of(" \n")   , 
                    boost::token_compress_on );

      S << size_of_elements + element_counter;

      // node number
      out += S.str();
      out += " ";
      S.str("");S.clear();

      // element property
      out += words[1];
      out += " ";

      // unstretched length
      out += words[2];
      out += " ";

      // anchor node
      lower = boost::lexical_cast<int>( words[3] ) + I*(size_of_nodes); // size_of_elements+1
      S << lower;
      out += S.str();
      out += " " ;
      S.str("");S.clear();

      // fairlead node
      upper = boost::lexical_cast<int>( words[4] ) + I*(size_of_nodes);
      S << upper;
      out += S.str();
      out += " " ;
      S.str("");S.clear();
	    
      // node flags. j=5 at start since the element option flag begin at words[5].
      for( unsigned int j=5 ; j<words.size() ; j++ ){
        out += words[j];
        out += " ";
      };

      out += "\n";
      element_data.push_back( out );      
      out.clear();
	    
      element_counter++;
    }
  } else {
    // ERROR! An exception cannot yet be thrown in this version of MAP
    // @todo: throw something here
  }
};


// ============================================================================
// GetCableLibraryData
//
// ============================================================================
std::string &MAP_InitInputType_class::
GetCableLibraryData( const unsigned int i )
{
  if ( i >= cable_library_data.size( ) ) {
    throw MAP_ERROR_4;
  };// END if
  
  return cable_library_data[i]; 
};


// ============================================================================
// GetNodeData
//
// ============================================================================
std::string &MAP_InitInputType_class::
GetNodeData( const unsigned int i ) 
{ 
  if ( i >= node_data.size( ) ) {
    throw MAP_ERROR_14;
  };// END if
  
  return node_data[i];          
};


// ============================================================================
// GetElementData
//
// ============================================================================
std::string &MAP_InitInputType_class::
GetSolverOptions( int i ) 
{ 
  return solver_options[i];       
};


// ============================================================================
// GetElementData 
//
// ============================================================================
std::string &MAP_InitInputType_class::
GetElementData( int i ) 
{ 
  return element_data[i];       
};


// ============================================================================
// SizeofCableLibrary
//
// ============================================================================
int MAP_InitInputType_class::
SizeofCableLibrary( ) const 
{ 
  return cable_library_data.size(); 
};


// ============================================================================
// SizeofSolverOptions
//
// ============================================================================
int MAP_InitInputType_class::
SizeofSolverOptions( ) const
{ 
  return solver_options.size(); 
};


// ============================================================================
// SizeofNodeData
//
// ============================================================================
int MAP_InitInputType_class::
SizeofNodeData( ) const
{ 
  return node_data.size();          
};


// ============================================================================
// SizeofElementData
//
// ============================================================================
int MAP_InitInputType_class::
SizeofElementData( ) const 
{ 
  return element_data.size();       
};


// ============================================================================
// SetDepth
//
// ============================================================================
void MAP_InitInputType_class::
SetDepth( const std::string &text )
{
  depth = text;
};


// ============================================================================
// SetGravity
//
// ============================================================================
void MAP_InitInputType_class::
SetGravity( const std::string &text ) 
{
  gravity = text;
};


// ============================================================================
// SetSeaDensity
//
// ============================================================================
void MAP_InitInputType_class::
SetSeaDensity( const std::string &text )
{
  sea_density = text;
};


// ============================================================================
// GetCableLibraryData
//
// ============================================================================
std::string &MAP_InitInputType_class::GetDepth( ) { 
  return depth; 
};


// ============================================================================
// GetCableLibraryData
//
// ============================================================================
std::string &MAP_InitInputType_class::
GetGravity( )
{  
  return gravity; 
};


// ============================================================================
// GetCableLibraryData
//
// ============================================================================
std::string &MAP_InitInputType_class::
GetSeaDensity( ) 
{ 
  return sea_density; 
};


// ============================================================================
// SetDepth
//
// Overloaded function for Fortran binding
// ============================================================================
void MAP_InitInputType_class::
SetDepth( const double text ) 
{
  this->depth_dbl = text;
};


// ============================================================================
// SetGravity
//
// ============================================================================
void MAP_InitInputType_class::
SetGravity( const double text ) 
{
  this->gravity_dbl = text;
};


// ============================================================================
// SetSeaDensity
//
// ============================================================================
void MAP_InitInputType_class::
SetSeaDensity( const double text )
{
  this->sea_density_dbl = text;
};


// ============================================================================
// GetDepth_double
//
// ============================================================================
double MAP_InitInputType_class::
GetDepth_double( ) const 
{ 
  return depth_dbl;     
}


// ============================================================================
// GetGravity_double
//
// ============================================================================
double MAP_InitInputType_class::
GetGravity_double( ) const 
{ 
  return gravity_dbl;   
};


// ============================================================================
// GetSeaDensity_double
//
// ============================================================================
double MAP_InitInputType_class::
GetSeaDensity_double( ) const 
{ 
  return sea_density_dbl; 
};
