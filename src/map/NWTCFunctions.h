/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   NWTCFunctions.h
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
