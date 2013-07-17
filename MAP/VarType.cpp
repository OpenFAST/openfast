/**
 * ====================================================================================================
 *                              VarType.cpp
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


#include "VarType.h"


/**
 * ====================================================================================================
 * setVarType 
 *
 * Static member
 * 
 * This functions check for the following:
 *    -- is the variable begins with a '#' in the MAP input file
 *    -- handles exceptions
 *    -- converts string arguments to doubles
 *
 * @input : T            -- VarType corresponding to FX, FY. FZ, M, B, X, Y, Z or Lu
 * @input : input_string --
 * @input : Error        -- Error numerical code
 * @input : Msg          -- Error message (string) that will eventually be passed back to the calling module 
 * ====================================================================================================
 */
void VarType::
setGenericVarType ( VarType           &T            , 
                    const std::string &input_string ,
                    MAP_ErrStat       &Error        , 
                    MAP_Message       &Msg) 
{
  std::string error_output = "";
    
  try {
    // if varaible starts with '#' :
    if (input_string[0] == '#') {
      // temp is a the equaivalent to 'in' but with '#' removed
      std::string temp = "";
      
      // initial guess flag for MAP
      if ( boost::iequals(input_string , "#") ) {
        temp = "9999.9";
      }
      else {
        unsigned int i = 1;
      
        // Now cycle through the string and write it to the This parameter,
        // but ignore the firest character (which will be a '#')
        while(i < input_string.size() ){
          temp += input_string[i];
          i++;
        };// END while
      };//END if
            
      T.value = boost::lexical_cast<double>( temp );	
      T.is_fixed = false;
    }
    else { // if variable does not start with a '#' :
      // we do not need to change the state of is_fixed
      // since it is initialized to true by default
      T.value = boost::lexical_cast<double>( input_string );
    };// END if
  } catch ( boost::bad_lexical_cast const& ) {
    std::string str = "";
    str += boost::lexical_cast < std::string > ( MAP_ERROR_26 );
    std::ostringstream S;
    S << T.index+1;
    str += "] : " + 
      MAP_ERROR_CODE_to_string.at( MAP_ERROR_26 ) + 
      T.name + 
      "in Node or Element " + 
      S.str() + 
      ".";

    Error.set_error_key( MAP_ERROR );

    // Here is the exception message
    Msg.RecordErrorToErrorList( str );

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
 * ====================================================================================================
 * writeGenericVarType 
 *
 * Prints the VarType varaible.name and VarType.value to the screen and profile.map file in the following format:
 *  
 *             'X[1] :  (-9999.9)'
 *
 * @input : VarType -- VarType variable we 
 *
 * @output : VarType.name + VarType.value in string form
 *
 * @todo  : Fix assert statements which are currently commented out
 * @todo  : out.size() < 10 is arbitrarily chosen. Fix this so the string can exceed a length of 10 
 * ====================================================================================================
 */
std::string VarType::
writeGenericVarType( VarType &T )
{
  // initialize the string we are writting for the outputs
  std::string out = "";
  std::ostringstream S;
  S << std::fixed << std::setprecision(3);

  // If T is an iterated value we are solving for, do this :
  if( T.is_fixed == false) {	
    if ( T.value == 9999.9 ) S << std::fixed << std::setprecision(1);

    // create the string we are writting as an output
    out += "   ";
    out += T.name;
    out += "[";
    S << T.index;
    out += S.str();
    S.str("");S.clear();

    // add extra spaces if the T.name varable is less than two
    // characters in length
    if ( T.name=="FX" || T.name=="FY" || T.name=="FZ" ) {
      out +="]:";
    }
    else {
      out +="]: ";
    };// END if

    do { // fill string with white spaces to bring it to the correct length
      out += " ";
    } while( out.size()<10 ); // 10 is arbitrarily chosen. 
                         
    S << T.value;
    out += "(";
    out += S.str();
    out += ")";//_TEXT_COLOR::END;
    S.str("");S.clear();

    //to be sure the length of 'out' is not greater than
    //length we are allow to across, this assert statement 
    //is raised. This only because visiable if the condition
    //is violated
    //	
    //@todo : fix this assert statement and restore it so it functions
    //        as usable code
    // assert( (out.size()-_TEXT_COLOR::STR_CHR)<_TEXT_COLOR::STR_LEN );

    // Fill the 'out' variable with white space if it is 
    // shorter than the space we want the variable to 
    // occupy
    do {
      out += " ";
    } while( (out.size() )<_TEXT_COLOR::STR_LEN);
  }
  else { // If the parameter is a fixed value, do this :
    // create the string we are writting as an output
    out += "   ";
    out += T.name;
    out += "[";
    S << T.index;
    out += S.str();
    S.str("");S.clear();

    // add extra spaces if the T.name varable is less than two
    // characters in length
    if ( T.name=="FX" || T.name=="FY" || T.name=="FZ" ) {
      out +="]:";
    }
    else {
      out +="]: ";
    };// END if

    do {
      out += " ";
    } while( out.size()<10 ); // @todo: 10 is arbitrarily chosen.                                  
    // Find a more permanent method
 
    S << T.value;
    out += S.str();
    S.str("");S.clear();
    
    // to be sure the length of 'out' is not greater than
    // length we are allow to across, this assert statement 
    // is raised. This only because visiable if the condition
    // is violated
    //	assert( out.size()<_TEXT_COLOR::STR_LEN );  
        
    // Fill the 'out' variable with white space if it is 
    // shorter than the space we want the variable to 
    // occupy
    do {
      out += " ";
    } while( out.size()<_TEXT_COLOR::STR_LEN );
  }//END if

  return out;
};


/**
 * ====================================================================================================
 * writeNodeData 
 *
 * Prints the VarType variable.name and VarType.value to the screen and profile.map file in the following format:
 *  
 *             '(-9999.9)'
 *
 * @input : VarType -- VarType variable we 
 *
 * @output : VarType.name + VarType.value in string form
 *
 * @todo  : Fix assert statements which are currently commented out
 * @todo  : out.size() < 10 is arbitrarily chosen. Fix this so the string can exceed a length of 10 
 * ====================================================================================================
 */
std::string VarType::
writeNodeData( VarType &T )
{
  // initlize the string we are writting for the outputs
  std::string out = "";
  std::ostringstream S;
  S << std::fixed << std::setprecision(3);

  // If T is an iterated value we are solving for, do this :
  if( T.is_fixed == false) {	
    if ( T.value == 9999.9 ) S << std::fixed << std::setprecision(1);
	                         
    S << T.value;
    out += "(";
    out += S.str();
    out += ")";//_TEXT_COLOR::END;
    S.str("");S.clear();
	
    // to be sure the length of 'out' is not greater than
    // length we are allow to across, this assert statement 
    // is raised. This only because visiable if the condition
    // is violated

    // Fill the 'out' variable with white space if it is 
    // shorter than the space we want the variable to 
    // occupy
    do {
      out += " ";
    } while( (out.size() )<_TEXT_COLOR::STR_LEN);
  }
    
  // If the parameter is a fixed value, do this :
  else {
    S << T.value;
    out += S.str();
    S.str("");S.clear();
	    
    // to be sure the length of 'out' is not greater than
    // length we are allow to across, this assert statement 
    // is raised. This only because visiable if the condition
    // is violated

    // Fill the 'out' variable with white space if it is 
    // shorter than the space we want the variable to 
    // occupy
    do {
      out += " ";
    } while( out.size()<_TEXT_COLOR::STR_LEN );
  }//END if

  return out;
};




/**
 * ====================================================================================================
 * writeGenericVarType 
 *
 * Prints the VarType varaible.name and VarType.value to the screen and profile.map file in 
 * the following format:
 *  
 *             'X[1] :           '
 *
 * @input : in           --
 * @input : T            -- VarType corresponding to FX, FY. FZ, M, B, X, Y, Z or Lu
 *
 * @output : string      -- VarType.name + VarType.value in string form
 *
 * @todo  : Fix assert statements which are currently commented out
 * @todo  : out.size() < 15 is arbitrarily chosen. Fix this so the string can exceed a length of 15 
 * ====================================================================================================
 */
std::string VarType::
writeGenericVarType_name( int i      , 
                          VarType &T )
{
  // initlize the string we are writting for the outputs
  std::string out = "";
  std::ostringstream S;
  S << std::fixed << std::setprecision(3);

  // create the string we are writting as an output
  out += T.name;
  out += "[";
  S << i+1;
  out += S.str();
  out += "]";
  S.str("");S.clear();
	
  // make sure the string array is long enough to put at least one white space 
  // at the tail; otherwise, the do-while loop will result in a program crash
  assert( out.size() < 14 );

  // fill string with white spaces to bring it to the correct length
  do {
    out += " ";
  } while( out.size()<15 ); // @note : 15 is arbitrarily chosen. 
 
  return out;
};


/**
 * ====================================================================================================
 * writeGenericVarType_value
 *
 * prints the VarType.value parameter in this style: 'X[1] :           '
 *
 * @input : i            -- the node/element index number
 * @input : T            -- VarType corresponding to FX, FY. FZ, M, B, X, Y, Z or Lu
 *
 * @output : string      -- the VarType.value as a std::string type 
 *
 * @todo  : out.size() < 15 is arbitrarily chosen. Fix this so the string can exceed a length of 15 
 * ====================================================================================================
 */
std::string VarType::
writeGenericVarType_value( int i      , 
                           VarType &T )
{
  // initlize the string we are writting for the outputs
  std::string out = "";
  std::ostringstream S;
  S << std::fixed << std::setprecision(3);

  // create the string we are writting as an output
  S << T.value;
  out += S.str();
  S.str("");S.clear();
	
  // make sure the string array is long enough to put at least one white space 
  // at the tail; otherwise, the do-while loop will result in a program crash
  assert( out.size() < 14 );

  // fill string with white spaces to bring it to the correct length
  do {
    out += " ";
  } while( out.size()<15 ); // @note : 15 is arbitrarily chosen. 
 
  return out;
};


/**
 * ====================================================================================================
 * SetUniversalErrorStat
 *
 * Records the error to the MAP_Message contect, and records and error status
 * ====================================================================================================
 */
static void 
SetUniversalErrorStat( MAP_ERROR_CODE code ,
                       MAP_ErrStat &Error  ,
                       MAP_Message &Msg  )
{
  std::string str = "";                                 
  str += boost::lexical_cast < std::string > ( code );  
  if ( str.size( ) == 1 ) str = " " + str;              
  str += "] : " + MAP_ERROR_CODE_to_string.at( code );  
      
  Error.set_error_key( MAP_ERROR );                     
  Msg.RecordErrorToErrorList( str );                                
      
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
