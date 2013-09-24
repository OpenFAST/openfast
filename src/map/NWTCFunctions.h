/**
 * ====================================================================================================
 *                              NWTCFunctions.h
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


#ifndef _NWTC_FUNCTIONS_H
#define _NWTC_FUNCTIONS_H


                                                          
void MSQS_Init(MAP_InitInputType_class       &Init     ,   
               MAP_InputType_class           &In       ,   
               MAP_ParameterType_class       &Param    ,   
               void* NULL_ContinuousState              ,         
               void* NULL_DiscreteState                ,         
               MAP_ConstraintStateType_class &Constrnt ,   
               MAP_OtherStateType_class      &Data     ,   
               MAP_OutputType_class          &Out      ,   
               void* NULL_TimeInterval                 ,         
               MAP_InitOutputType_class      &InitOut  ,
               MAP_ErrStat_class             &Error    ,   
               MAP_Message_class             &Msg      );       

                                                          
void MSQS_UpdateStates( float                         T        ,
                        int                           interval ,
                        MAP_InputType_class           &In      ,   
                        MAP_ParameterType_class       &Param   ,   
                        void* NULL_ContinuousState             ,        
                        void* NULL_DiscreteState               ,        
                        MAP_ConstraintStateType_class &Constrnt,  
                        MAP_OtherStateType_class      &Data    ,  
                        MAP_ErrStat_class             &Error   ,  
                        MAP_Message_class             &Msg     );     

                                                                  
void MSQS_CalcOutput( float                         T         ,          
                      MAP_InputType_class           &In       ,    
                      MAP_ParameterType_class       &Param    ,    
                      void* NULL_ContinuousState              ,          
                      void* NULL_DiscreteState                ,          
                      MAP_ConstraintStateType_class &Constrnt ,    
                      MAP_OtherStateType_class      &Data     ,    
                      MAP_OutputType_class          &Out      ,    
                      MAP_ErrStat_class             &Error    ,    
                      MAP_Message_class             &Msg      );       

                                                                  
void MSQS_End( MAP_InputType_class           &In       ,           
               MAP_ParameterType_class       &Param    ,           
               void* NULL_ContinuousState              ,                 
               void* NULL_DiscreteState                ,                 
               MAP_ConstraintStateType_class &Constrnt ,           
               MAP_OtherStateType_class      &Data     ,           
               MAP_OutputType_class          &Out      ,           
               MAP_ErrStat_class             &Error    ,           
               MAP_Message_class             &Msg      );              
                                                             

void output_python( MAP_OtherStateType_class &Data , 
                    MAP_Message_class        &Msg  ); 



#endif // _NWTC_FUNCTIONS_H
