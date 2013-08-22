/**
 * ====================================================================================================
 *                              Catenary.h
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
