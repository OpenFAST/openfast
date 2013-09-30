/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   CatenaryEquation.h
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


#ifndef _CATENARYEQUATIONS_H
#define _CATENARYEQUATIONS_H


class Element;  // class predefintion. Class is declared below.

/**
 * ==========   Solving functions   ===============     <--------------------+
 * Catenary equations                                            //          |
 *                                                               //          |
 * For each element, there are two equations we                  //          |
 * are solving for:                                              //          |
 *  -- Vertical catenary equation.                               //          |
 *  -- Horizontal catenary equation.                             //          |
 *                                                               //          |
 * These equations are defined below as nested classes           //          |
 */                                                              //          |
                                                                 //          |
struct HORIZONTAL_CATENARY_EQ {                                  //          |
    Element * const this_element; // can't change reference      //          |
                                  // to this_element             //          |
    double Fh;      // = this_element->H.value;                  //          |
    double Fv;      // = this_element->V.value;                  //          |
    double omega;   // = this_element->omega;                    //          |
    double Lu;      // = this_element->Lu.value;                 //          |
    double length;  // = this_element->h;                        //          |
    double EA;      // = this_element->GetEA()
    double Lb;      // = length laying on sea bed                //          |
    double Cb;      // = cable friction coefficient              //          |
    double mu;      // = cable mass per length                   //          |
    double operator()();                                         //          |
HORIZONTAL_CATENARY_EQ( Element &aP ) : this_element( &aP ) {}   //          |
};                                                               //          |
                                                                 //          |
                                                                 //          |
struct VERTICAL_CATENARY_EQ {                                    //          |
    Element * const this_element; // can't change reference      //          |
                                  // to this_element             //          |
    double Fh;      // = this_element->H.value;                  //          |
    double Fv;      // = this_element->V.value;                  //          |
    double omega;   // = this_element->omega;                    //          |
    double Lu;      // = this_element->Lu.value;                 //          |
    double height;  // = this_element->h;                        //          |
    double EA;      // = this_element->GetEA()
    double Lb;      // = length laying on sea bed                //          |
    double Cb;      // = cable friction coefficient              //          |
    double mu;      // = cable mass per length                   //          |
    double operator()();                                         //          |
VERTICAL_CATENARY_EQ( Element &aP ) : this_element( &aP ) {}     //          |
};                                                               //  --------+
//============================================================================


/**
 * ==========   Catenary equation derivatives   ======     <-----------------+
 *                                                               //          |  
 */                                                              //          |                                                           
                                                                 //          |
struct DXDH {                                                    //          |
    Element * const this_element; // can't change reference      //          |
                                  // to this_element             //          |
    double Fh;      // = this_element->H.value;                  //          |
    double Fv;      // = this_element->V.value;                  //          |
    double omega;   // = this_element->omega;                    //          |
    double Lu;      // = this_element->Lu.value;                 //          |
    double height;  // = this_element->h;                        //          |
    double EA;
    double Lb;      // = length laying on sea bed                //          |
    double Cb;      // = cable friction coefficient              //          |
    double mu;      // = cable mass per length                   //          |
    double operator()();                                         //          |
DXDH( Element &aP ) : this_element( &aP ) {}                     //          |
};                                                               //          |
                                                                 //          |
                                                                 //          |
struct DXDV {                                                    //          |
    Element * const this_element; // can't change reference      //          |
                                  // to this_element             //          |
    double Fh;      // = this_element->H.value;                  //          |
    double Fv;      // = this_element->V.value;                  //          |
    double omega;   // = this_element->omega;                    //          |
    double Lu;      // = this_element->Lu.value;                 //          |
    double height;  // = this_element->h;                        //          |
    double EA;
    double Lb;      // = length laying on sea bed                //          |
    double Cb;      // = cable friction coefficient              //          |
    double mu;      // = cable mass per length                   //          |
    double operator()();                                         //          |
DXDV( Element &aP ) : this_element( &aP ) {}                     //          |
};                                                               //          |
                                                                 //          |
                                                                 //          |
struct DZDH {                                                    //          |
    Element * const this_element; // can't change reference      //          |
                                  // to this_element             //          |
    double Fh;      // = this_element->H.value;                  //          |
    double Fv;      // = this_element->V.value;                  //          |
    double omega;   // = this_element->omega;                    //          |
    double Lu;      // = this_element->Lu.value;                 //          |
    double height;  // = this_element->h;                        //          |
    double EA;
    double Lb;      // = length laying on sea bed                //          |
    double Cb;      // = cable friction coefficient              //          |
    double mu;      // = cable mass per length                   //          |
    double operator()();                                         //          |
DZDH( Element &aP ) : this_element( &aP ) {}                     //          |
};                                                               //          |
                                                                 //          |
                                                                 //          |
struct DZDV {                                                    //          |
    Element * const this_element; // can't change reference      //          |
                                  // to this_element             //          |
    double Fh;      // = this_element->H.value;                  //          |
    double Fv;      // = this_element->V.value;                  //          |
    double omega;   // = this_element->omega;                    //          |
    double Lu;      // = this_element->Lu.value;                 //          |
    double height;  // = this_element->h;                        //          |
    double EA;
    double Lb;      // = length laying on sea bed                //          |
    double Cb;      // = cable friction coefficient              //          |
    double mu;      // = cable mass per length                   //          |
    double operator()();                                         //          |
DZDV( Element &aP ) : this_element( &aP ) {}                     //          |
};                                                               //  --------+
//============================================================================


#endif // _CATENARYEQUATIONS_H
