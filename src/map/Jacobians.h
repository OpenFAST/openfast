/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   Jacobians.h
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
 * 
 * Base class for the A, B and -B^T Jacobian blocks in J.
 * 
 * @see       UserData::SetJacB( )
 * @see       UserData::GetJacA( )
 * @property  xf       X fairlead position (global)
 * @property  xa       X anchor position (global)
 * @property  yf       Y fairlead position (global)
 * @property  ya       Y anchor position (global)
 * @property  sign     -1 or +1
 * @property  F        H.value (horizontal force)
 * @property  Element  point to element (presverse states without passing arguments)
 * @default   return 0.0
 * 
 */                                                        
                                                           
struct A_DERIVS {                                          
  Element *this_element; // can't change reference         
                                                           
  double xf;    // x fairlead displacement                    
  double xa;    // x anchor displacement                      
  double yf;    // y fairlead displacement                    
  double ya;    // y anchor displacement                      
  double sign;  // + or -, depending on cable attachment.
  double F;     // Horizontal force
                                                           
A_DERIVS( ) : xf   ( 0.0 ) ,                               
    xa             ( 0.0 ) ,                               
    yf             ( 0.0 ) ,                               
    ya             ( 0.0 ) ,                               
    sign           ( 1.0 ) ,                               
    F              ( 0.0 ) {}                              
                                                           
  void setElement( Element *T )      { this_element = T;  }
  void setSign   ( const double in ) { this->sign   = in; }
                                                           
  virtual double getDeriv() { return 0.0; }                
};                                                         
                                                           
/**                                                        
 *             d Hx                                        
 * dHxdxf     ------                                       
 *             d xf                                        
 */                                                        
struct dHxdxf : A_DERIVS{                                  
  double getDeriv() {                                      
    xf = this_element->GetXf();                            
    xa = this_element->GetXa();                            
    yf = this_element->GetYf();                            
    ya = this_element->GetYa();                            
    F  = this_element->GetH();                             
                                                           
//    return sign*( F )*( (yf-ya)*(yf-ya) ) / pow( ( (xf-xa)*(xf-xa) + (yf-ya)*(yf-ya) ) , 1.5 );   
    return sign*( F )*pow( (yf-ya) , 2) / pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );   
  }                                                      
};                                                       
                                                         
                                                         
/**                                                      
 *             d Hx                                      
 * dHxdxa     ------                                     
 *             d xa                                      
 */                                                      
struct dHxdxa : A_DERIVS{                                
  double getDeriv() {                                    
    xf = this_element->GetXf();                          
    xa = this_element->GetXa();                          
    yf = this_element->GetYf();                          
    ya = this_element->GetYa();                          
    F  = this_element->GetH();                           
                                                         
//    return sign*( F )*( (yf-ya)*(yf-ya) ) / pow( ( (xf-xa)*(xf-xa) + (yf-ya)*(yf-ya) ) , 1.5 );   
    return sign*( F )*pow( (yf-ya) , 2) / pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );   
  }                                                      
};                                                       
                                                         
/**                                                      
 *             d Hx                                      
 * dHxdyf     ------                                     
 *             d yd                                      
 */                                                      
struct dHxdyf : A_DERIVS{                                
  double getDeriv() {                                    
    xf = this_element->GetXf();                          
    xa = this_element->GetXa();                          
    yf = this_element->GetYf();                          
    ya = this_element->GetYa();                          
    F  = this_element->GetH();                           
	                                                 
//    return -sign*( F )*( (xf-xa)*(yf-ya) ) / pow( ( (xf-xa)*(xf-xa) + (yf-ya)*(yf-ya) ) , 1.5 ); 
    return -sign*( F )*( (xf-xa)*(yf-ya) ) / pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );   
  }                                                      
};                                                       
                                                         
                                                         
/**                                                      
 *             d Hx                                      
 * dHxdya     ------                                     
 *             d ya                                      
 */                                                      
struct dHxdya : A_DERIVS{                                
  double getDeriv() {                                    
    xf = this_element->GetXf();                          
    xa = this_element->GetXa();                          
    yf = this_element->GetYf();                          
    ya = this_element->GetYa();                          
    F  = this_element->GetH();                           

//    return -sign*( F )*( (xf-xa)*(yf-ya) ) / pow( ( (xf-xa)*(xf-xa) + (yf-ya)*(yf-ya) ) , 1.5 );
    return -sign*( F )*( (xf-xa)*(yf-ya) ) / pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );   
  }                                                      
};                                                       
                                                         
                                                         
/**                                                      
 *             d Hy                                      
 * dHydxf     ------                                     
 *             d xf                                      
 */                                                      
struct dHydxf : A_DERIVS{                                
  double getDeriv() {                                    
    xf = this_element->GetXf();                          
    xa = this_element->GetXa();                          
    yf = this_element->GetYf();                          
    ya = this_element->GetYa();                          
    F  = this_element->GetH();                           
    
//    return -sign*( F )*( (xf-xa)*(yf-ya) ) / pow( ( (xf-xa)*(xf-xa) + (yf-ya)*(yf-ya) ) , 1.5 ); 
    return -sign*( F )*( (xf-xa)*(yf-ya) ) / pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );   
  }                                                      
};                                                       
                                                         
                                                         
/**                                                      
 *             d Hy                                      
 * dHydxa     ------                                     
 *             d xa                                      
 */                                                      
struct dHydxa : A_DERIVS{                                
  double getDeriv() {                                    
    xf = this_element->GetXf();                          
    xa = this_element->GetXa();                          
    yf = this_element->GetYf();                          
    ya = this_element->GetYa();                          
    F  = this_element->GetH();                           
                                                         
//    return -sign*( F )*( (yf-ya)*(xf-xa) ) / pow( ( (xf-xa)*(xf-xa) + (yf-ya)*(yf-ya) ) , 1.5 );
    return -sign*( F )*( (yf-ya)*(xf-xa) ) / pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );   
  }                                                      
};                                                       
                                                         
                                                         
/**                                                      
 *             d Hy                                      
 * dHydyf     ------                                     
 *             d yf                                      
 */                                                      
struct dHydyf : A_DERIVS{                                
  double getDeriv() {                                    
    xf = this_element->GetXf();                          
    xa = this_element->GetXa();                          
    yf = this_element->GetYf();                          
    ya = this_element->GetYa();                          
    F  = this_element->GetH();                           

    // error here
//    return sign*( F )*( (xf-xa)*(xf-xa) ) / pow( ( (xf-xa)*(xf-xa) + (yf-ya)*(yf-ya) ) , 1.5 );  
    return sign*( F )*pow( (xf-xa) , 2 ) / pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );  
  }                                                      
};                                                       
                                                         
                                                         
/**                                                      
 *             d Hy                                      
 * dHydya     ------                                     
 *             d ya                                      
 */                                                      
struct dHydya : A_DERIVS{                                
  double getDeriv() {                                    
    xf = this_element->GetXf();                          
    xa = this_element->GetXa();                          
    yf = this_element->GetYf();                          
    ya = this_element->GetYa();                          
    F  = this_element->GetH();                           

//    return sign*( F )*( (xf-xa)*(xf-xa) ) / pow( ( (xf-xa)*(xf-xa) + (yf-ya)*(yf-ya) ) , 1.5 );
    return sign*( F )*pow( (xf-xa) , 2 ) / pow( ( pow(xf-xa,2) + pow(yf-ya,2) ) , 1.5 );     
  }                                                        
};                                                         
//============================================================================


/**
 * A_BLOCK
 */
struct A_BLOCK {
private:
  int ix;
  int iy;
  std::vector <A_DERIVS_ptr> Acell;          
    
public:    

  void initializeACell( Element *elem_in , CASE_DERIV in , const double polarity );
  void setRow( const int i ) { this->iy = i; }
  void setCol( const int i ) { this->ix = i; }
  int row( ) const { return  iy; }
  int col( ) const { return  ix; }  
  
  double GetDerivativeForBlockA() { 
    double sum = 0.0;	
    for( unsigned int i=0 ; i<Acell.size() ; i++ )  sum += Acell[i]->getDeriv();
    return sum;
  }
};


/**
 * B_BLOCK
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
