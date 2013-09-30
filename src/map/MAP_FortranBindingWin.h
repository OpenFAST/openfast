/**
 * ====================================================================================================
 *                              FortranBinding.h
 * ====================================================================================================
 *	     
 * Copyright Sept. 2012
 * 
 * Author: Marco Masciola 
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

#ifdef _WIN64
  #include <windows.h>
  #define CALL __declspec( dllexport ) //define something for Windows (64-bit)
#elif _WIN32
  #include <windows.h>
  #define CALL __declspec( dllexport ) //define something for Windows (32-bit)
#elif __APPLE__
  // Unsupported platform
#elif __linux
  #define CALL __attribute__((dllexport) )
#elif __unix // all unices not caught above
  #define CALL __attribute__((dllexport) )
#elif __posix
    // POSIX
#endif


void 
MAPSetInputArray( MAP_InputType            *input , 
                  MAP_InputType_class      *u     ,
                  MAP_OtherStateType       *other ,
                  MAP_OtherStateType_class *o );

void 
MAPSetParameterArray( MAP_ParameterType        *param , 
                      MAP_ParameterType_class  *p     ,
                      MAP_OtherStateType       *other ,
                      MAP_OtherStateType_class *o );

void 
MAPSetConstraintArray( MAP_ConstraintStateType       *constr , 
                       MAP_ConstraintStateType_class *z      ,
                       MAP_OtherStateType            *other  ,
                       MAP_OtherStateType_class      *o );

void 
MAPSetOtherArray( void*                    NULL_A,
                  void*                    NULL_B,
                  MAP_OtherStateType       *other  ,
                  MAP_OtherStateType_class *o );

void 
MAPSetOutputArray( MAP_OutputType           *output , 
                   MAP_OutputType_class     *y      ,
                   MAP_OtherStateType       *other  ,
                   MAP_OtherStateType_class *o );

void 
MAPUnpackConstraint( MAP_ConstraintStateType       *constrData ,
                     MAP_ConstraintStateType_class *z          ,
                     MAP_OtherStateType            *otherState );

void 
MAPUnpackInput( MAP_InputType       *inputData ,
                MAP_InputType_class *u         ,
                MAP_OtherStateType  *otherState );

void 
MAPUnpackOther( MAP_OtherStateType       *otherState ,
                MAP_OtherStateType_class *o);

void 
MAPUnpackParameter( MAP_ParameterType       *paramData ,
                    MAP_ParameterType_class *p         ,
                    MAP_OtherStateType      *otherState );
void 
MAPPackOutput( MAP_OutputType       *output ,
               MAP_OutputType_class *y      ,
               MAP_OtherStateType   *otherState );

void 
MAPPackParameter( MAP_ParameterType       *paramData ,
                  MAP_ParameterType_class *p         ,
                  MAP_OtherStateType      *otherState );

void 
MAPPackConstraint( MAP_ConstraintStateType       *constrData ,
                   MAP_ConstraintStateType_class *z         ,
                   MAP_OtherStateType             *otherState );
void 
MAPPackInput( MAP_InputType       *input ,
              MAP_InputType_class *u     ,
              MAP_OtherStateType  *otherState );


template <class Type> void
MAPPack( Type   *ptr     ,
         int    *index   ,
         int    len      ,
         double *array   ,
         int    *counter )
{
  for(int k=0 ; k<len ; k++ ){                      
    array[k] = ptr->GetVar( index[*counter] );                   
    (*counter)++;                                                    
  };
};

/**
 * ==========   SetArrayIndex   ===============
 *                                                            
 * These function set the indexing number for the Input, Outputs, 
 * Constraints, ect. vectors that are inside MAP. To change a variables 
 * in MAP, Object->SetVar( index , value ) is called. These function 
 * set index. This is called by the FAST glue code 
 */                                                           
template <class Type> void 
MAPSetArrayIndex( const std::string   &argStr , 
                  int                 &len    , 
                  Type                *Ptr    ,                       
                  std::vector<int>    &index  , 
                  std::vector<double> &value )
{                           
  len = 0;
  value.clear();
  for ( int i=0 ; i<Ptr->size() ; i++ ) {
    if ( argStr == Ptr->GetElementName(i) ) {
      index.push_back(i);
      value.push_back( Ptr->GetVar(i) );
      len++;
    };
  };
};//  ==========   END SetArrayIndex   ===============


template <class Type> void
MAPUnpack( Type   *ptr     ,
           int    *index   ,
           int    len      ,
           double *array   ,
           int    *counter )
{
  for(int k=0 ; k<len ; k++ ){         
    ptr->SetVar( index[*counter], array[k] );
    (*counter)++;
  }                                            
};

/**
 * ==========   Delete Arrays in Types   ===============     
 *
 * Used in the CALLMAP_MSQS_End routine to clean up allocated 
 * arrays. These functions are called once.                   
 */                                                           
void                                                          
MAPDeleteOutputArray( MAP_OutputType *OutData )               
{                                                             
  delete[] OutData->FX;                                       
  delete[] OutData->FY;                                       
  delete[] OutData->FZ;                                       
};                                                            
                                                              

void                                                          
MAPDeleteInputArray( MAP_InputType *InputData )               
{                                                             
  delete[] InputData->X;                                      
  delete[] InputData->Y;                                      
  delete[] InputData->Z;                                      
};                                                            
                                                              

void                                                          
MAPDeleteConstraintArray( MAP_ConstraintStateType *ConstrData )
{
  delete[] ConstrData->X;                                    
  delete[] ConstrData->Y;                                    
  delete[] ConstrData->Z;                                    
  delete[] ConstrData->H;                                    
  delete[] ConstrData->V;                                    
};                                                           
                                                             

void                                                         
MAPDeleteParameterArray( MAP_ParameterType *ParamData )      
{                                                            
  delete[] ParamData->Diam;                                  
  delete[] ParamData->MassDenInAir;                          
  delete[] ParamData->EA;                                    
  delete[] ParamData->CB;                                    
  delete[] ParamData->Lu;                                    
  delete[] ParamData->X;                                     
  delete[] ParamData->Y;                                     
  delete[] ParamData->Z;                                     
  delete[] ParamData->FX;                                    
  delete[] ParamData->FY;                                    
  delete[] ParamData->FZ;                                    
  delete[] ParamData->M;                                     
  delete[] ParamData->B;                                     
};                                                           
                                                             
void                                                         
MAPDeleteOtherArray( MAP_OtherStateType *OtherData )         
{                                                            
  delete[] OtherData->FX;                                    
  delete[] OtherData->FY;                                    
  delete[] OtherData->FZ;                                    
  delete[] OtherData->o_index;                               
  delete[] OtherData->FLAGS_index;                         
  delete[] OtherData->PLOT_flag;                           
  delete[] OtherData->X_POS_flag;                          
  delete[] OtherData->Y_POS_flag;                          
  delete[] OtherData->Z_POS_flag;                          
  delete[] OtherData->X_FORCE_flag;                        
  delete[] OtherData->Y_FORCE_flag;                        
  delete[] OtherData->Z_FORCE_flag;                        
  delete[] OtherData->LINE_TENSION_flag;                   
  delete[] OtherData->OMIT_CONTACT_flag;                   
  delete[] OtherData->LAY_LENGTH_flag;                     
                                                             
                                                             
  // inputs                                                  
  delete[] OtherData->u_index;                               
                                                             
  // Output                                                  
  delete[] OtherData->y_index;                               
                                                             
  // Constraint                                              
  delete[] OtherData->z_index;                               
                                                             
  // Parameters                                              
  delete[] OtherData->p_index;                               
                                                             
  // Other States                                            
};                                                           
//=================================================================================================
