/**
 * ====================================================================================================
 *                              UserData.cpp
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


#include "UserData.h"


/**
 * ====================================================================================================
 * sizeOfConstraint
 * 
 * @output :
 * ====================================================================================================
 */
int UserData::sizeOfConstraint ( ) const {
  return this->constraint->size();
};


/**
 * ====================================================================================================
 * setConstraint
 *
 *  
 * ====================================================================================================
 */
void UserData::setConstraint( int i , double p ) {
  this->constraint->setVar( i , p);
};


/**
 * ====================================================================================================
 * getConstraint 
 *
 * @output :
 * ====================================================================================================
 */
double UserData::getConstraint( const int i ) const {
  return this->constraint->getVar( i );
};


/**
 * ====================================================================================================
 * sizeOfNode
 * 
 * @output : size of the nodes identified as 'Connect'. We are bglecting any 'Vessel' or 'Fix'
 *           nodes since those do not need a Newton equation.
 * ====================================================================================================
 */
int UserData::sizeOfNode( ) { 
  return this->node.size( ); 
};


/**
 * ====================================================================================================
 * sizeOfElement
 *
 * @output : number of elements. This is use dto identify the number of equations being solved for the
 *           continuous analytical cable equations.
 * ====================================================================================================
 */
int UserData::sizeOfElement( ) { 
  return this->element.size(); 
};


/**
 * ====================================================================================================
 * getNumElemEqs
 * 
 * @output :
 * ====================================================================================================
 */
int UserData::getNumElemEqs( ) {
  int cnt = 0;
  for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {
    cnt++;
    cnt++; 
  };//END for
  return cnt;
};


/**
 * ====================================================================================================
 * getNumNodeEqs
 *
 * count number of element equations
 * 
 * @note : This is hard coded assuming there are two equations per element.
 *         This is not necessarily true, as in the case of a vertical cable.
 * 
 * @todo : Find a solution to the above problem.
 * 
 * @output :
 * ====================================================================================================
 */
int UserData::getNumNodeEqs( ) {
  int cnt = 0;
  for ( int i=0 ; i<this->sizeOfNode() ; i++ ) {	
    if (this->node[i]->getXNewtonEquationFlag()==true) cnt++;
    if (this->node[i]->getYNewtonEquationFlag()==true) cnt++;
    if (this->node[i]->getZNewtonEquationFlag()==true) cnt++;
  };//END for	
  return cnt;
};


/**
 * ====================================================================================================
 * setJacB
 * ====================================================================================================
 */
void UserData::setJacB( int indexX , int indexY , int elem_index , DERIV dd , double polarity ){

  B_BLOCK_ptr b_block_ptr( new B_BLOCK );      // create new memory on the heap (dynamically allocated)

  b_block_ptr->setElementReference( this->element[ elem_index ] );
  b_block_ptr->setIx              ( indexX                      );
  b_block_ptr->setIy              ( indexY                      );
  b_block_ptr->setPartial         ( dd                          );
  b_block_ptr->setSign            ( polarity                    );

  JacB.push_back( b_block_ptr );
};


/**
 * ====================================================================================================
 * setDeriv
 * ====================================================================================================
 */
void A_BLOCK::initializeACell( Element *elem_in , CASE_DERIV in , const double polarity ) { 

  A_DERIVS_ptr a_derivs_ptr( new A_DERIVS );  // create new memory on the heap (dynamically allocated)
    
  switch( in ){
  case DHXDXF :
    a_derivs_ptr.reset( new dHxdxf );
    break;
  case DHXDXA :
    a_derivs_ptr.reset( new dHxdxa );
    break;
  case DHXDYF :
    a_derivs_ptr.reset( new dHxdyf );
    break;
  case DHXDYA : 
    a_derivs_ptr.reset( new dHxdya );
    break;
  case DHYDXF : 
    a_derivs_ptr.reset( new dHydxf );
    break;
  case DHYDXA : 
    a_derivs_ptr.reset( new dHydxa );
    break;
  case DHYDYF :
    a_derivs_ptr.reset( new dHydyf ); 
    break;
  case DHYDYA : 
    a_derivs_ptr.reset( new dHydya );
    break;
  };// END case

  a_derivs_ptr->setElement( elem_in );
  a_derivs_ptr->setSign( polarity );

  Acell.push_back( a_derivs_ptr );
  //std::cout << "Size of Acell : " << Acell[0]->getDeriv() << std::endl;
};


/**
 * ====================================================================================================
 * setJacA
 * ====================================================================================================
 */
void UserData::setJacA( const int row , const int i , const DERIV dd ) {

  A_BLOCK_ptr a_block_ptr( new A_BLOCK );  // create new memory on the heap (dynamically allocated)

  for ( int j=0 ; j<this->sizeOfElement() ; j++ ) {                                      
    if( this->element[j]->compareNodeAddressWithFairlead( *node[i] ) ) {
      switch( dd ) {
      case DX :
        a_block_ptr->initializeACell( this->element[j] , DHXDXF , 1.0 );
        a_block_ptr->setRow( row );
        a_block_ptr->setCol( row );
        break;
      case DY :
        a_block_ptr->initializeACell( this->element[j] , DHYDYF , 1.0 );
        a_block_ptr->setRow( row );
        a_block_ptr->setCol( row );
        break;
      case DZ : 
        a_block_ptr->setRow( row );
        a_block_ptr->setCol( row );
        break;
      default :
        std::cout << "Error in createDeriv" << std::endl;
      };// END switch
    }

    else if( this->element[j]->compareNodeAddressWithAnchor( *node[i] ) ) {  
      switch( dd ) {
      case DX :
        a_block_ptr->initializeACell( this->element[j] , DHXDXA , 1.0 );
        a_block_ptr->setRow( row );
        a_block_ptr->setCol( row );		
        break;
      case DY :
        a_block_ptr->initializeACell( this->element[j] , DHYDYA , 1.0 );
        a_block_ptr->setRow( row );
        a_block_ptr->setCol( row );
        break;
      case DZ : 		
        a_block_ptr->setRow( row );
        a_block_ptr->setCol( row );
        break;
      default :
        std::cout << "Error in createDeriv" << std::endl;
      };// END switch
    };// END if-else
  };// END for

  JacA.push_back( a_block_ptr );        

  switch( dd ) {
  case DX :
    if ( this->node[i]->getYNewtonEquationFlag()==true ) {
      A_BLOCK_ptr b_block_ptr( new A_BLOCK );
      for ( int j=0 ; j<this->sizeOfElement() ; j++ ) {                                      
        if( this->element[j]->compareNodeAddressWithFairlead( *node[i] ) ) {
          b_block_ptr->initializeACell( this->element[j] , DHXDYF , 1.0 );
          b_block_ptr->setRow( row+1 );
          b_block_ptr->setCol( row   );
        }

        else if( this->element[j]->compareNodeAddressWithAnchor( *node[i] ) ) {
          b_block_ptr->initializeACell( this->element[j] , DHXDYA , 1.0 );
          b_block_ptr->setRow( row+1 );
          b_block_ptr->setCol( row   );
        }
      }
      JacA.push_back( b_block_ptr );        
    }
    break;

  case DY :
    if ( this->node[i]->getXNewtonEquationFlag()==true ) {
      A_BLOCK_ptr b_block_ptr( new A_BLOCK );
      for ( int j=0 ; j<this->sizeOfElement() ; j++ ) {                                      
        if( this->element[j]->compareNodeAddressWithFairlead( *node[i] ) ) {
          b_block_ptr->initializeACell( this->element[j] , DHYDXF , 1.0 );
          b_block_ptr->setRow( row-1 );
          b_block_ptr->setCol( row   );
        }

        else if( this->element[j]->compareNodeAddressWithAnchor( *node[i] ) ) {
          b_block_ptr->initializeACell( this->element[j] , DHYDXA , 1.0 );
          b_block_ptr->setRow( row-1 );
          b_block_ptr->setCol( row   );
        }
      }
      JacA.push_back( b_block_ptr );        
    }
    break;

  case DZ :	
    break;

  default :
    std::cout << "Error in createDeriv" << std::endl;
  };// END switch
};







void UserData::findIndex( const int nodeIndex , const int col , const DERIV dd ){
    
  int cnt=0;
  for ( int j=0 ; j<this->sizeOfNode() ; j++ ) {					
    if ( this->node[j]->getXNewtonEquationFlag()==true ) {  		    
      for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {                   
        if( this->element[i]->compareNodeAddressWithFairlead( *node[nodeIndex] ) ) {
          if( this->element[i]->compareNodeAddressWithAnchor( *node[j] ) ) {
            A_BLOCK_ptr a_block_ptr( new A_BLOCK );  
			
            if(dd==DX) a_block_ptr->initializeACell( this->element[i] , DHXDXF , -1.0 );
            if(dd==DY) a_block_ptr->initializeACell( this->element[i] , DHYDXF , -1.0 );
			
            a_block_ptr->setRow( cnt );
            a_block_ptr->setCol( col );
            JacA.push_back( a_block_ptr );        
          }			
        }

        if( this->element[i]->compareNodeAddressWithAnchor( *node[nodeIndex] ) ) {
          if( this->element[i]->compareNodeAddressWithFairlead( *node[j] ) ) {
            A_BLOCK_ptr a_block_ptr( new A_BLOCK );  

            if(dd==DX) a_block_ptr->initializeACell( this->element[i] , DHXDXA , -1.0 );
            if(dd==DY) a_block_ptr->initializeACell( this->element[i] , DHYDXA , -1.0 );

            a_block_ptr->setRow( cnt );
            a_block_ptr->setCol( col );
            JacA.push_back( a_block_ptr );
          }			
        }

      }
      cnt++;
    }	    
	
    if ( this->node[j]->getYNewtonEquationFlag()==true ) {  		    
      for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {                   
        if( this->element[i]->compareNodeAddressWithFairlead( *node[nodeIndex] ) ) {
          if( this->element[i]->compareNodeAddressWithAnchor( *node[j] ) ) {
            A_BLOCK_ptr a_block_ptr( new A_BLOCK );  
                
            if(dd==DX) a_block_ptr->initializeACell( this->element[i] , DHXDYF , -1.0 );
            if(dd==DY) a_block_ptr->initializeACell( this->element[i] , DHYDYF , -1.0 );

            a_block_ptr->setRow( cnt );
            a_block_ptr->setCol( col );
            JacA.push_back( a_block_ptr );       
          }			
        }

        if( this->element[i]->compareNodeAddressWithAnchor( *node[nodeIndex] ) ) {
          if( this->element[i]->compareNodeAddressWithFairlead( *node[j] ) ) {
            A_BLOCK_ptr a_block_ptr( new A_BLOCK );  

            if(dd==DX) a_block_ptr->initializeACell( this->element[i] , DHXDYA , -1.0 );
            if(dd==DY) a_block_ptr->initializeACell( this->element[i] , DHYDYA , -1.0 );

            a_block_ptr->setRow( cnt );
            a_block_ptr->setCol( col );
            JacA.push_back( a_block_ptr );       
          }			
        }

      }
      cnt++;
    }	    

    if ( this->node[j]->getZNewtonEquationFlag()==true ) {  		    
      cnt++;
    }
  }
};


/**
 * ====================================================================================================
 * initializeJacobian
 *
 * Find the non-zero patter for:
 * 
 *   J = [  A     B ]
 *       [ -B^T   C ]
 *
 * Note that the C block is not defined here because it is trivial
 * ====================================================================================================
 */
void UserData::initializeJacobian(){
    
  int m=0;
  int n=0;

  /**
   * =======  Initialize A block  ======     <-------------------------------------------------------------+
   *                                                                                           //          |
   * We are finding the non-zero pattern for the B matrix bock in J.                           //          |
   * The address of each element node is compared with the address of the                      //          |
   * node solved in the newton equation.                                                       //          |
   */                                                                                          //          |
  for ( int i=0 ; i<this->sizeOfNode() ; i++ ) {	                                         //          |
    if ( this->node[i]->getXNewtonEquationFlag()==true ) {                                 //          |
      this->setJacA( m , i , DX );                
      m++;
    };                                                                

    if ( this->node[i]->getYNewtonEquationFlag()==true ) {  
      this->setJacA( m , i , DY );                
      m++;
    };                                                                

    if ( this->node[i]->getZNewtonEquationFlag()==true ) {  
      this->setJacA( m , i , DZ );                
      m++;
    };                                                                
  };                                                   




  int cnt=0;    
    
  for ( int k=0 ; k<this->sizeOfNode() ; k++ ) {
	
    if ( this->node[k]->getXNewtonEquationFlag()==true ) {  
      this->findIndex( k , cnt , DX );
      cnt++;
    }
    if ( this->node[k]->getYNewtonEquationFlag()==true ) {  
      this->findIndex( k , cnt , DY );
      cnt++;
    }
    if ( this->node[k]->getZNewtonEquationFlag()==true ) {  
      cnt++;
    }
	
  };// END for
    
    
    //============== <END> ===================================================================================

    
    /**
     * =======  Initialize B block  ======     <-------------------------------------------------------------+
     *                                                                                           //          |
     * We are finding the non-zero pattern for the B matrix bock in J.                           //          |
     * The address of each element node is compared with the address of the                      //          |
     * node solved in the newton equation.                                                       //          |
     */                                                                                          //          |
  m=0;                                                                                         //          | 
  for ( int i=0 ; i<this->sizeOfNode() ; i++ ) {	                                         //          |
    /**                                                                                      //          |
     * =======  Non-zero for (dFx)/(dH) and (dFx)/(dV)  ======     <-------------------+     //          |
     */                                                                          //    |     //          |
    //    |     //          |
    n=0;                                                                         //    |     //          |
    if (this->node[i]->getXNewtonEquationFlag()==true) {                         //    |     //          |
      for ( int j=0 ; j<this->sizeOfElement() ; j++ ) {                        //    |     //          |
        if( this->element[j]->compareNodeAddressWithFairlead( *node[i] ) ) { //    |     //          |
          this->setJacB( n , m , j , DX , -1.0 );                          //    |     //          |
        };                                                                   //    |     //          |
        //    |     //          |
        if( this->element[j]->compareNodeAddressWithAnchor( *node[i] ) ) {   //    |     //          |
          this->setJacB( n , m , j , DX , 1.0 );                           //    |     //          |
        };                                                                   //    |     //          |
        n+=2;                                                                //    |     //          |
      };                                                                       //    |     //          |
      m++;                                                                     //    |     //          |
    };                                                                           // ---+     //          |
    //============== <END> =============================================================     //          |
    //          |
    /**                                                                                      //          |
     * =======  Non-zero for (dFy)/(dH) and (dFy)/(dV)  ======     <-------------------+     //          |
     */                                                                          //    |     //          |
    //    |	 //          |
    n=0;                                                                         //    |     //          |
    // solve Y direction Newton equation                                         //    |     //          |
    if (this->node[i]->getYNewtonEquationFlag()==true) {                         //    |     //          |
      for ( int j=0 ; j<this->sizeOfElement() ; j++ ) {                        //    |     //          |
        if( this->element[j]->compareNodeAddressWithFairlead( *node[i] ) ) { //    |     //          |
          this->setJacB( n , m , j , DY , -1.0 );                          //    |     //          |
        };                                                                   //    |     //          |
        //    |     //          |
        if( this->element[j]->compareNodeAddressWithAnchor( *node[i] ) ) {   //    |     //          |
          this->setJacB( n , m , j , DY , 1.0 );                           //    |     //          |
        };                                                                   //    |     //          |
        n+=2;                                                                //    |     //          |
      };                                                                       //    |     //          |
      m++;                                                                     //    |     //          |
    };                                                                           // ---+     //          |
    //============== <END> =============================================================     //          |
    //          |
    //          |
    /**                                                                                      //          |
     * =======  Non-zero for (dFz)/(dH) and (dFz)/(dV)  ======     <-------------------+     //          |
     */                                                                          //    |     //          |
    //    |     //          |
    n=1;                                                                         //    |     //          |
    // solve Z direction Newton equation                                         //    |     //          |
    if (this->node[i]->getZNewtonEquationFlag()==true) {                         //    |     //          |
      for ( int j=0 ; j<this->sizeOfElement() ; j++ ) {                        //    |     //          |
        if( this->element[j]->compareNodeAddressWithFairlead( *node[i] ) ) { //    |     //          |
          this->setJacB( n , m , j , DZ , -1.0 );                          //    |     //          |
        };                                                                   //    |     //          |
        //    |     //          |
        if( this->element[j]->compareNodeAddressWithAnchor( *node[i] ) ) {   //    |     //          |
          this->setJacB( n , m , j , DZ , 1.0 );                           //    |     //          |
        };                                                                   //    |     //          |
        n+=2;                                                                //    |     //          |
      };                                                                       //    |     //          |
      m++;                                                                     //    |     //          |
    };                                                                           // ---+     //          |
    //============== <END> =============================================================     //          |
  };                                                                                           //   -------+
  //============== <END> ===================================================================================
};


/**
 * ====================================================================================================
 * B_BLOCK()
 * 
 * The B matrix clock is based purely on the geometry of how the lines are connected to one another
 *
 *  J = [  A     B ]
 *      [ -B^T   C ]
 * 
 * Get the derivatives used to populate the B block in the Jacobian.
 * ====================================================================================================
 */
double B_BLOCK::getDeriv(){
//double B_BLOCK::operator()() {
  if ( partial==DX ) {
    //return 0;
    return sign*cos( this->elem->getPsi() ); 
  }
  else if ( partial==DY ){
    //return 0;
    return sign*sin( this->elem->getPsi() ); 
  }
  else if ( partial==DZ ){
    //return 0;
    return sign;
  }
  else {
    std::cout << "MAP : Error in B_BLOCK() function object" << std::endl; 
    return 0.0;
  };
};
