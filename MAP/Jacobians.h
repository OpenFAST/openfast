/**
 * ====================================================================================================
 *                              Jacobians.h
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


#ifndef _JACOBIANS_H
#define _JACOBIANS_H


enum DERIV { NONE , DX , DY , DZ , FAIRLEAD , ANCHOR };
enum CASE_DERIV{ DHXDXF , DHXDXA , DHXDYF , DHXDYA , DHYDXF , DHYDXA , DHYDYF , DHYDYA };

struct A_BLOCK;
struct A_DERIVS;
struct B_BLOCK;

typedef boost::shared_ptr <A_BLOCK>  A_BLOCK_ptr;
typedef boost::shared_ptr <A_DERIVS> A_DERIVS_ptr;
typedef boost::shared_ptr <B_BLOCK>  B_BLOCK_ptr;


/**
 * ==========   A Jacobian Block base   ===============     <----------------+
 *                                                               //          |
 *                                                               //          |
 */                                                              //          |
                                                                 //          |
struct A_DERIVS {                                                //          |
  Element *this_element; // can't change reference               //          |
                                                                 //          |
  double xf; // x fairlead displacement                          //          |
  double xa; // x anchor displacement                            //          |
  double yf; // y fairlead displacement                          //          |
  double ya; // y anchor displacement                            //          |
  double sign;                                                   //          |
  double F;                                                      //          |
                                                                 //          |
A_DERIVS( ) : xf   ( 0.0 ) ,                                     //          |
    xa             ( 0.0 ) ,                                     //          |
    yf             ( 0.0 ) ,                                     //          |
    ya             ( 0.0 ) ,                                     //          |
    sign           ( 1.0 ) ,                                     //          |
    F              ( 0.0 ) {}                                    //          |
                                                                 //          |
  void setElement( Element *T )      { this_element = T;  }      //          |
  void setSign   ( const double in ) { this->sign   = in; }      //          |
                                                                 //          |
  virtual double getDeriv() { return 0.0; }                      //          |
};                                                               //          |
                                                                 //          |
/**                                                              //          |
 * --------------------------------------------------            //          |
 *             d Hx                                              //          |
 * dHxdxf     ------                                             //          |
 *             d xf                                              //          |
 */                                                              //          |
struct dHxdxf : A_DERIVS{                                        //          |
  double getDeriv() {                                            //          |
    xf = this_element->getXf();                                  //          |
    xa = this_element->getXa();                                  //          |
    yf = this_element->getYf();                                  //          |
    ya = this_element->getYa();                                  //          |
    F  = this_element->getH();                                   //          |
                                                                 //          |
    return sign*( F )*pow( (yf-ya) , 2)                          //          |
      /pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );             //          |
  }                                                              //          |
};                                                               //          |
                                                                 //          |
                                                                 //          |
/**                                                              //          |
 * --------------------------------------------------            //          |
 *             d Hx                                              //          |
 * dHxdxa     ------                                             //          |
 *             d xa                                              //          |
 */                                                              //          |
struct dHxdxa : A_DERIVS{                                        //          |
  double getDeriv() {                                            //          |
    xf = this_element->getXf();                                  //          |
    xa = this_element->getXa();                                  //          |
    yf = this_element->getYf();                                  //          |
    ya = this_element->getYa();                                  //          |
    F  = this_element->getH();                                   //          |
                                                                 //          |
    return sign*( F )*pow( (yf-ya) , 2)                          //          |
      /pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );             //          |
  }                                                              //          |
};                                                               //          |
                                                                 //          |
/**                                                              //          |
 * --------------------------------------------------            //          |
 *             d Hx                                              //          |
 * dHxdyf     ------                                             //          |
 *             d yd                                              //          |
 */                                                              //          |
struct dHxdyf : A_DERIVS{                                        //          |
  double getDeriv() {                                            //          |
    xf = this_element->getXf();                                  //          |
    xa = this_element->getXa();                                  //          |
    yf = this_element->getYf();                                  //          |
    ya = this_element->getYa();                                  //          |
    F  = this_element->getH();                                   //          |
	                                                         //          |
    return -sign*( F )*( (xf-xa)*(yf-ya) )                       //          |
      /pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );             //          |
  }                                                              //          |
};                                                               //          |
                                                                 //          |
                                                                 //          |
/**                                                              //          |
 * --------------------------------------------------            //          |
 *             d Hx                                              //          |
 * dHxdya     ------                                             //          |
 *             d ya                                              //          |
 */                                                              //          |
struct dHxdya : A_DERIVS{                                        //          |
  double getDeriv() {                                            //          |
    xf = this_element->getXf();                                  //          |
    xa = this_element->getXa();                                  //          |
    yf = this_element->getYf();                                  //          |
    ya = this_element->getYa();                                  //          |
    F  = this_element->getH();                                   //          |
                                                                 //          |
    return -sign*( F )*( (xf-xa)*(yf-ya) )                       //          |
      /pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );             //          |
  }                                                              //          |
};                                                               //          |
                                                                 //          |
                                                                 //          |
/**                                                              //          |
 * --------------------------------------------------            //          |
 *             d Hy                                              //          |
 * dHydxf     ------                                             //          |
 *             d xf                                              //          |
 */                                                              //          |
struct dHydxf : A_DERIVS{                                        //          |
  double getDeriv() {                                            //          |
    xf = this_element->getXf();                                  //          |
    xa = this_element->getXa();                                  //          |
    yf = this_element->getYf();                                  //          |
    ya = this_element->getYa();                                  //          |
    F  = this_element->getH();                                   //          |
                                                                 //          |
    return -sign*( F )*( (xf-xa)*(yf-ya) )                       //          |
      /pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );             //          |
  }                                                              //          |
};                                                               //          |
                                                                 //          |
                                                                 //          |
/**                                                              //          |
 * --------------------------------------------------            //          |
 *             d Hy                                              //          |
 * dHydxa     ------                                             //          |
 *             d xa                                              //          |
 */                                                              //          |
struct dHydxa : A_DERIVS{                                        //          |
  double getDeriv() {                                            //          |
    xf = this_element->getXf();                                  //          |
    xa = this_element->getXa();                                  //          |
    yf = this_element->getYf();                                  //          |
    ya = this_element->getYa();                                  //          |
    F  = this_element->getH();                                   //          |
                                                                 //          |
    return -sign*( F )*( (yf-ya)*(xf-xa) )                       //          |
      /pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );             //          |
  }                                                              //          |
};                                                               //          |
                                                                 //          |
                                                                 //          |
/**                                                              //          |
 * --------------------------------------------------            //          |
 *             d Hy                                              //          |
 * dHydyf     ------                                             //          |
 *             d yf                                              //          |
 */                                                              //          |
struct dHydyf : A_DERIVS{                                        //          |
  double getDeriv() {                                            //          |
    xf = this_element->getXf();                                  //          |
    xa = this_element->getXa();                                  //          |
    yf = this_element->getYf();                                  //          |
    ya = this_element->getYa();                                  //          |
    F  = this_element->getH();                                   //          |
                                                                 //          |
    return sign*( F )*pow( (xf-xa) , 2 )                         //          |
      /pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );             //          |
  }                                                              //          |
};                                                               //          |
                                                                 //          |
                                                                 //          |
/**                                                              //          |
 * --------------------------------------------------            //          |
 *             d Hy                                              //          |
 * dHydya     ------                                             //          |
 *             d ya                                              //          |
 */                                                              //          |
struct dHydya : A_DERIVS{                                        //          |
  double getDeriv() {                                            //          |
    xf = this_element->getXf();                                  //          |
    xa = this_element->getXa();                                  //          |
    yf = this_element->getYf();                                  //          |
    ya = this_element->getYa();                                  //          |
    F  = this_element->getH();                                   //          |
                                                                 //          |
    return sign*( F )*pow( (xf-xa) , 2 )                         //          |
      /pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );             //          |
  }                                                              //          |
};                                                               //  --------+
//============================================================================


/**
 * ====================================================================================================
 * A_BLOCK
 *
 * ====================================================================================================
 */
struct A_BLOCK {
private:
  int ix;
  int iy;
  std::vector <A_DERIVS_ptr> Acell;          
    
public:    

  //void createDeriv( const int i , const int j , Element* elem_in , DERIV xyz , DERIV location );
    
  void setRow( const int i ) { this->iy = i; }
  void setCol( const int i ) { this->ix = i; }

  int row( ) const { return  iy; }
  int col( ) const { return  ix; }

  void   initializeACell( Element *elem_in , CASE_DERIV in , const double polarity );
    
  double getDeriv() { 
    double sum = 0.0;	
    for( unsigned int i=0 ; i<Acell.size() ; i++ )  sum += Acell[i]->getDeriv();
    return sum;
  }
};

/**
 * ====================================================================================================
 * B_BLOCK
 *
 * ====================================================================================================
 */
struct B_BLOCK {
private:
  Element *elem;
  DERIV partial;
  int ix;
  int iy;    
  double sign;

public:

  void setElementReference( Element *elem_in ) { this->elem = elem_in; }
  void setIx              ( const int i      ) { this->ix = i;         }
  void setIy              ( const int i      ) { this->iy = i;         }
  void setPartial         ( const DERIV in   ) { this->partial = in;   }
  void setSign            ( const double i   ) { this->sign = i;       }

  int row() const { return  iy; }
  int col() const { return  ix; }
  double getDeriv();
};



#endif // _JACOBIANS_H
