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
 * This is a static member, meaning it can be called throughout MAP just by calling
 * VarType::SetGenericVarType(...)
 * 
 * This functions check for the following:
 *   - is the variable begins with a '#' in the MAP input file
 *   - handles exceptions
 *   - converts string arguments to doubles
 *
 * @param : var           VarType pointer corresponding to FX, FY. FZ, M, B, X, Y, Z or Lu
 * @param : input_string  string to set the name of the variable
 * @param : err           Error enum code
 * @param : msg           Error message (string) that will eventually be passed back to 
 *                        the calling module 
 */
void VarType::
SetGenericVarType ( VarType           &var          , 
                    const std::string &input_string ,
                    MAP_ErrStat_class       &err          , 
                    MAP_Message_class       &msg          ) 
{
  std::string error_output = "";
  double tempVal = 0.0;
  std::string temp = "";
  
  try {    
    if (input_string[0] == '#') {  // if varaible starts with '#' :
      // initial guess flag for MAP      
      if ( boost::iequals(input_string , "#") ) {
        temp = "9999.9";
      } else {
        unsigned int i = 1;
        
        // Now cycle through the string and write it to the This parameter,
        // but ignore the firest character (which will be a '#')
        while(i < input_string.size() ){
          temp += input_string[i];
          i++;
        };
      };
  
      tempVal = boost::lexical_cast<double>( temp );	
      if ( fabs(tempVal)<=1e-4 ) {
        var.value = 0.0;
      } else {
        var.value = tempVal; 
      }
      var.is_fixed = false;
// @change mdm : this was removed to account for the case where the value is very close to zero. 
//               This was creating error per bug report
//      var.value = boost::lexical_cast<double>( temp );	
//      var.is_fixed = false;
// @endchange
    } else { // if variable does not start with a '#' :      
      // we do not need to change the state of is_fixed since it is initialized to true by default
      tempVal = boost::lexical_cast<double>( input_string );	
      if ( fabs(tempVal)<=1e-4 ) {
        var.value = 0.0;
      } else {
        var.value = tempVal; 
      }
// @change mdm : this was removed to account for the case where the value is very close to zero. 
//               This was creating error per bug report
//      var.value = boost::lexical_cast<double>( input_string );       
// @endchange
    }
  } catch ( boost::bad_lexical_cast const& ) {
    std::string str = "";
    std::ostringstream S;

    S << var.index+1;
    str += var.name + " in Node or Element " + S.str() + ".";
    
    MAPSetUniversalErrorStat( MAP_ERROR_26, str, err, msg );
  }
};


/**
 * Prints the VarType var.name and VarType var.value to the screen and profile.map file in 
 * the following format:
 *  
 *             'X[1] :  (-9999.9)'
 *
 * @param   var          VarType variable we writting 
 * @return  std::string  VarType.name + VarType.value in string form
 * @todo                 Fix assert statements which are currently commented out
 * @todo                 out.size() < 10 is arbitrarily chosen. Fix this so the 
 *                       string can exceed a length of 10 
 */
std::string VarType::
WriteGenericVarType( VarType &var )
{
  // initialize the string we are writting for the outputs
  std::string out = "";
  std::ostringstream S;
  S << std::fixed << std::setprecision(3);

  // If T is an iterated value we are solving for, do this :
  if( var.is_fixed == false) {	
    if ( var.value == 9999.9 ) S << std::fixed << std::setprecision(1);

    // create the string we are writting as an output
    S << var.index;
    out += "   " + var.name + "[" + S.str() + out += "]:";
    S.str("");S.clear();

    // add extra spaces if the T.name varable is less than two characters in length
    if ( var.name!="FX" || var.name!="FY" || var.name!="FZ" ) out +=" ";

    // fill string with white spaces to bring it to the correct length
    do { out += " "; } while( out.size()<10 ); // 10 is arbitrarily chosen. 
                         
    S << var.value;
    out += "(" +  S.str() + ")"; // adding "()" around value to indicate it is iterated
    S.str("");S.clear();

    // Fill the 'out' variable with white space if it is 
    // shorter than the space we want the variable to 
    // occupy
    do { out += " "; } while( (out.size() )<_TEXT_COLOR::STR_LEN);
  } else { // If the parameter is a fixed value, do this :
    // create the string we are writting as an output
    S << var.index;
    out += "   " + var.name + "[" + S.str() + "]:";

    S.str("");S.clear();

    // add extra spaces if the var.name varable is less than two
    // characters in length
    if ( var.name!="FX" || var.name!="FY" || var.name!="FZ" ) out +=" ";

    do { out += " ";} while( out.size()<10 ); // @todo: 10 is arbitrarily chosen.
 
    S << var.value;
    out += S.str();
    S.str("");S.clear();
            
    // Fill the 'out' variable with white space if it is 
    // shorter than the space we want the variable to 
    // occupy
    do { out += " "; } while( out.size()<_TEXT_COLOR::STR_LEN );
  }//END if

  return out;
};


// ====================================================================================================
// WriteNodeData 
//
// Prints the VarType variable.name and VarType.value to the screen and profile.map file in 
// the following format:
//  
//             '(-9999.9)'
//
// @input : VarType -- VarType variable we 
//
// @output : VarType.name + VarType.value in string form
//
// @todo  : Fix assert statements which are currently commented out
// @todo  : out.size() < 10 is arbitrarily chosen. Fix this so the string can exceed a length of 10 
// ====================================================================================================
std::string VarType::
WriteNodeData( VarType &var )
{
  // initlize the string we are writting for the outputs
  std::string out = "";
  std::ostringstream S;
  S << std::fixed << std::setprecision(3);

  // If var is an iterated value we are solving for, do this :
  if( var.is_fixed == false) {	
    if ( var.value == 9999.9 ) S << std::fixed << std::setprecision(1);
	                         
    S << var.value;
    out += "(" +  S.str() + ")";
    S.str("");S.clear();
	
    // Fill the 'out' variable with white space if it is 
    // shorter than the space we want the variable to occupy
    do { out += " "; } while( (out.size() )<_TEXT_COLOR::STR_LEN);
  } else { // If the parameter is a fixed value, do this :
    S << var.value;
    out += S.str();
    S.str("");S.clear();

    // Fill the 'out' variable with white space if it is 
    // shorter than the space we want the variable to 
    // occupy
    do { out += " "; } while( out.size()<_TEXT_COLOR::STR_LEN );
  }//END if

  return out;
};


// ====================================================================================================
// WriteGenericVarType 
//
// Prints the VarType varaible.name and VarType.value to the screen and profile.map file in 
// the following format:
//  
//             'X[1] :           '
//
// @input : in           --
// @input : var          -- VarType corresponding to FX, FY. FZ, M, B, X, Y, Z or Lu
//
// @output : string      -- VarType.name + VarType.value in string form
//
// @todo  : Fix assert statements which are currently commented out
// @todo  : out.size() < 15 is arbitrarily chosen. Fix this so the string can exceed a length of 15 
// ====================================================================================================
std::string VarType::
WriteGenericVarType_name( const int i  , 
                          VarType &var )
{
  // initlize the string we are writting for the outputs
  std::string out = "";
  std::ostringstream S;
  S << std::fixed << std::setprecision(3);

  // create the string we are writting as an output
  S << i+1;
  out += var.name + "[" + S.str() + "]";
  S.str("");S.clear();
	
  // make sure the string array is long enough to put at least one white space 
  // at the tail; otherwise, the do-while loop will result in a program crash
  assert( out.size() < 14 );

  // fill string with white spaces to bring it to the correct length
  do { out += " "; } while( out.size()<15 ); // @note : 15 is arbitrarily chosen. 
 
  return out;
};


// ====================================================================================================
// WriteGenericVarType_value
//
// prints the VarType.value parameter in this style: 'X[1] :           '
//
// @input : i            -- the node/element index number
// @input : var          -- VarType corresponding to FX, FY. FZ, M, B, X, Y, Z or Lu
//
// @output : string      -- the VarType.value as a std::string type 
//
// @todo  : out.size() < 15 is arbitrarily chosen. Fix this so the string can exceed a length of 15 
// ====================================================================================================
std::string VarType::
WriteGenericVarType_value( const int i  , 
                           VarType &var )
{
  // initlize the string we are writting for the outputs
  std::string out = "";
  std::ostringstream S;
  S << std::fixed << std::setprecision(3);

  // create the string we are writting as an output
  S << var.value;
  out += S.str();
  S.str("");S.clear();
	
  // make sure the string array is long enough to put at least one white space 
  // at the tail; otherwise, the do-while loop will result in a program crash
  assert( out.size() < 14 );

  // fill string with white spaces to bring it to the correct length
  do { out += " "; } while( out.size()<15 ); // @note : 15 is arbitrarily chosen. 
 
  return out;
};


// ====================================================================================================
// MAPSetUniversalErrorStat
//
// Records the error to the MAP_Message_class contect, and records and error status
// ====================================================================================================
void 
MAPSetUniversalErrorStat( MAP_ERROR_CODE code     ,
                          std::string    &userStr ,
                          MAP_ErrStat_class    &err     ,
                          MAP_Message_class    &msg     )
{
  std::string str = "";                                 
  std::ostringstream S;                                 

  str += boost::lexical_cast < std::string > ( code );  
  str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + userStr;  

  if ( code >= MAP_WARNING_1 ) {
    err.set_error_key( MAP_WARNING ); // MAP did not quite fail
    msg.RecordToWarningList( str );   // Let users know what the error is
  } else {
    err.set_error_key( MAP_ERROR ); // MAP failed
    msg.RecordToErrorList( str );   // Let users know why MAP failed
  }
      
  str.erase(2,1);                                       
  S.str("");
  S.clear();                                  
  
  S << ">>>> " + str + "    :    MAP error in file ";     
  S << __FILE__;                                        
  S << " on line ";                                     
  S << __LINE__;                                        
  S << " occurred.";                                    
  
  msg.WriteErrorToOutputFile( S.str() );                           
}
