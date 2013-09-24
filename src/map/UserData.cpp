/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                        UserData.cpp
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


#include "UserData.h"


// ====================================================================================================
// sizeOfConstraint
// 
// @output :
// ====================================================================================================
int UserData::
sizeOfConstraint ( ) const 
{
  return this->constraint->size();
};


// ====================================================================================================
// setConstraint
//
//  
// ====================================================================================================
void UserData::
setConstraint( const int    index , 
               const double value ) 
{
  this->constraint->SetVar( index, value );
};


// ====================================================================================================
// getConstraint 
//
// @output :
// ====================================================================================================
double UserData::
getConstraint( const int index ) const 
{
  return this->constraint->GetVar( index );
};


// ====================================================================================================
// sizeOfNode
// 
// @output : size of the nodes identified as 'Connect'. We are bglecting any 'Vessel' or 'Fix'
//           nodes since those do not need a Newton equation.
// ====================================================================================================
int UserData::
sizeOfNode( ) 
{ 
  return this->node.size( ); 
};


// ====================================================================================================
// sizeOfElement
//
// @output : number of elements. This is use dto identify the number of equations being solved for the
//           continuous analytical cable equations.
// ====================================================================================================
int UserData::
sizeOfElement( ) 
{ 
  return this->element.size(); 
};


// ====================================================================================================
// getNumElemEqs
// 
// @output :
// ====================================================================================================
int UserData::
getNumElemEqs( ) 
{
  int cnt = 0;
  for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {
    cnt++;
    cnt++; 
  }
  return cnt;
};


// ====================================================================================================
// getNumNodeEqs
//
// count number of element equations
// 
// @note : This is hard coded assuming there are two equations per element.
//         This is not necessarily true, as in the case of a vertical cable.
// 
// @todo : Find a solution to the above problem.
// 
// @output :
// ====================================================================================================
int UserData::
getNumNodeEqs( ) 
{
  int cnt = 0;
  for ( int i=0 ; i<this->sizeOfNode() ; i++ ) {	
    if (this->node[i]->GetXNewtonEquationFlag()==true) cnt++;
    if (this->node[i]->GetYNewtonEquationFlag()==true) cnt++;
    if (this->node[i]->GetZNewtonEquationFlag()==true) cnt++;
  }
  return cnt;
};


// ====================================================================================================
// setJacB
// ====================================================================================================
void UserData::
setJacB( const int    index_x    , 
         const int    index_y    , 
         const int    elem_index , 
         const DERIV  dd         , 
         const double polarity   )
{
  B_BLOCK_ptr b_block_ptr( new B_BLOCK );      // create new memory on the heap (dynamically allocated)

  b_block_ptr->setElementReference( this->element[ elem_index ] );
  b_block_ptr->setIx              ( index_x                     );
  b_block_ptr->setIy              ( index_y                    );
  b_block_ptr->setPartial         ( dd                          );
  b_block_ptr->setSign            ( polarity                    );

  JacB.push_back( b_block_ptr );
};


// ====================================================================================================
// setDeriv
// ====================================================================================================
void A_BLOCK::
initializeACell( Element      *elem_in , 
                 CASE_DERIV   in       , 
                 const double polarity ) 
{ 
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
  a_derivs_ptr->setSign( polarity );   // either '+' or '-'

  Acell.push_back( a_derivs_ptr );
};


// ====================================================================================================
// setJacA
// ====================================================================================================
void UserData::
setJacA( const int row   , 
         const int index , 
         const DERIV Di  ) 
{
  A_BLOCK_ptr a_block_ptr( new A_BLOCK );  // create new memory on the heap (dynamically allocated)

  /*
      Diagonal components of the A Jacobian block      <-------------------------------------------+

      Di controls the denominator in dFi/Di = {dFx/dx ; dFz/dz ; dFz/dz }
  */
  for ( int j=0 ; j<this->sizeOfElement() ; j++ ) {                                      
    if( this->element[j]->CompareNodeAddressWithFairlead( *node[index] ) ) { // if the node is a fairlead
      switch( Di ) {
      case DX : // dF/dX
        a_block_ptr->initializeACell( this->element[j] , DHXDXF , 1.0 );
        a_block_ptr->setRow( row ); // same rows and column since diagonal
        a_block_ptr->setCol( row );
        break;
      case DY : // dF/dY
        a_block_ptr->initializeACell( this->element[j] , DHYDYF , 1.0 );
        a_block_ptr->setRow( row ); // same rows and column since diagonal
        a_block_ptr->setCol( row );
        break;
      case DZ : // dF/dZ
        a_block_ptr->setRow( row ); // same rows and column since diagonal
        a_block_ptr->setCol( row );
        break;
      default :
        std::string str = "";
        MAPSetUniversalErrorStat( MAP_ERROR_82 , str, *errPtr, *msgPtr );    
      }
    } else if( this->element[j]->CompareNodeAddressWithAnchor( *node[index] ) ) { // if the node is an anchor
      switch( Di ) {
      case DX : // dF/dX
        a_block_ptr->initializeACell( this->element[j] , DHXDXA , 1.0 );
        a_block_ptr->setRow( row );
        a_block_ptr->setCol( row );		
        break;
      case DY : // dF/dY
        a_block_ptr->initializeACell( this->element[j] , DHYDYA , 1.0 );
        a_block_ptr->setRow( row );
        a_block_ptr->setCol( row );
        break;
      case DZ : // dF/dZ	
        a_block_ptr->setRow( row );
        a_block_ptr->setCol( row );
        break;
      default :
        std::string str = "";
        MAPSetUniversalErrorStat( MAP_ERROR_82 , str, *errPtr, *msgPtr );    
      }
    }
  }
  JacA.push_back( a_block_ptr );        

  /*
      Off-Diagonal components of the A Jacobian block    <---------------------------------------+
  */
  switch( Di ) {
  case DX :
    if ( this->node[index]->GetYNewtonEquationFlag()==true ) {
      A_BLOCK_ptr a_off_diag_block_ptr( new A_BLOCK );
      for ( int j=0 ; j<this->sizeOfElement() ; j++ ) {                                      
        if( this->element[j]->CompareNodeAddressWithFairlead( *node[index] ) ) {
          a_off_diag_block_ptr->initializeACell( this->element[j] , DHXDYF , 1.0 ); // this creates a new A_DERIV
          a_off_diag_block_ptr->setRow( row+1 );
          a_off_diag_block_ptr->setCol( row   );
        } else if( this->element[j]->CompareNodeAddressWithAnchor( *node[index] ) ) { // this creates a new A_DERIV
          a_off_diag_block_ptr->initializeACell( this->element[j] , DHXDYA , 1.0 );
          a_off_diag_block_ptr->setRow( row+1 );
          a_off_diag_block_ptr->setCol( row   );
        }
      }
      JacA.push_back( a_off_diag_block_ptr );        
    }
    break;
  case DY :
    if ( this->node[index]->GetXNewtonEquationFlag()==true ) {
      A_BLOCK_ptr a_off_diag_block_ptr( new A_BLOCK );
      for ( int j=0 ; j<this->sizeOfElement() ; j++ ) {                                      
        if( this->element[j]->CompareNodeAddressWithFairlead( *node[index] ) ) {
          a_off_diag_block_ptr->initializeACell( this->element[j] , DHYDXF , 1.0 );
          a_off_diag_block_ptr->setRow( row-1 );
          a_off_diag_block_ptr->setCol( row   );
        } else if( this->element[j]->CompareNodeAddressWithAnchor( *node[index] ) ) {
          a_off_diag_block_ptr->initializeACell( this->element[j] , DHYDXA , 1.0 );
          a_off_diag_block_ptr->setRow( row-1 );
          a_off_diag_block_ptr->setCol( row   );
        }
      }
      JacA.push_back( a_off_diag_block_ptr );        
    }
    break;
  case DZ :	
    break;
  default :
    std::string str = "";
    MAPSetUniversalErrorStat( MAP_ERROR_91 , str, *errPtr, *msgPtr );    
  }
};


// ====================================================================================================
//
// ====================================================================================================
void UserData::findIndex( const int nodeIndex , const int col , const DERIV dd ){
    
  int cnt=0;
  for ( int j=0 ; j<this->sizeOfNode() ; j++ ) {					
    if ( this->node[j]->GetXNewtonEquationFlag()==true ) {  		    
      for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {                   
        if( this->element[i]->CompareNodeAddressWithFairlead( *node[nodeIndex] ) ) {
          if( this->element[i]->CompareNodeAddressWithAnchor( *node[j] ) ) {
            A_BLOCK_ptr a_block_ptr( new A_BLOCK );  
			
            if(dd==DX) a_block_ptr->initializeACell( this->element[i] , DHXDXF , -1.0 );
            if(dd==DY) a_block_ptr->initializeACell( this->element[i] , DHYDXF , -1.0 );
			
            a_block_ptr->setRow( cnt );
            a_block_ptr->setCol( col );
            JacA.push_back( a_block_ptr );        
          }			
        }

        if( this->element[i]->CompareNodeAddressWithAnchor( *node[nodeIndex] ) ) {
          if( this->element[i]->CompareNodeAddressWithFairlead( *node[j] ) ) {
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
	
    if ( this->node[j]->GetYNewtonEquationFlag()==true ) {  		    
      for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {                   
        if( this->element[i]->CompareNodeAddressWithFairlead( *node[nodeIndex] ) ) {
          if( this->element[i]->CompareNodeAddressWithAnchor( *node[j] ) ) {
            A_BLOCK_ptr a_block_ptr( new A_BLOCK );  
                
            if(dd==DX) a_block_ptr->initializeACell( this->element[i] , DHXDYF , -1.0 );
            if(dd==DY) a_block_ptr->initializeACell( this->element[i] , DHYDYF , -1.0 );

            a_block_ptr->setRow( cnt );
            a_block_ptr->setCol( col );
            JacA.push_back( a_block_ptr );       
          }			
        }

        if( this->element[i]->CompareNodeAddressWithAnchor( *node[nodeIndex] ) ) {
          if( this->element[i]->CompareNodeAddressWithFairlead( *node[j] ) ) {
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

    if ( this->node[j]->GetZNewtonEquationFlag()==true ) {  		    
      cnt++;
    }
  }
};


// ====================================================================================================
// initializeJacobian
//
// Find the non-zero patter for:
// 
//   J = [  A     B ]
//       [ -B^T   C ]
//
// Note that the C block is not defined here because it is trivial
// ====================================================================================================
void UserData::
initializeJacobian( )
{
  int m   = 0; // rows
  int n   = 0; // columns
  int cnt = 0;        


  /*
      Initialize A block       <-------------------------------------------------------------+
      
      We are finding the non-zero pattern for the B matrix bock in J. The address of each 
      element node is compared with the address of the node solved in the newton equation.
  */
  for ( int i=0 ; i<this->sizeOfNode() ; i++ ) {
    if ( this->node[i]->GetXNewtonEquationFlag()==true ) {
      this->setJacA( m , i , DX );                
      m++;
    }
    
    if ( this->node[i]->GetYNewtonEquationFlag()==true ) {  
      this->setJacA( m , i , DY );                
      m++;
    }
    
    if ( this->node[i]->GetZNewtonEquationFlag()==true ) {  
      this->setJacA( m , i , DZ );                
      m++;
    }
  }

  for ( int k=0 ; k<this->sizeOfNode() ; k++ ) {	
    if ( this->node[k]->GetXNewtonEquationFlag()==true ) {  
      this->findIndex( k , cnt , DX );
      cnt++;
    }
    if ( this->node[k]->GetYNewtonEquationFlag()==true ) {  
      this->findIndex( k , cnt , DY );
      cnt++;
    }
    if ( this->node[k]->GetZNewtonEquationFlag()==true ) {  
      cnt++;
    }	
  }

    
  /*
      Initialize B block       <-------------------------------------------------------------+
                                                                     
      We are finding the non-zero pattern for the B matrix bock in J.      
      The address of each element node is compared with the address of the 
      node solved in the newton equation.                                  
  */
  m = 0;                                                                    
  for ( int i=0 ; i<this->sizeOfNode() ; i++ ) {	                  


    /*
        Non-zero for (dFx)/(dH) and (dFx)/(dV)  

        Solve X direction Newton equation
    */
    n = 0;
    if ( this->node[i]->GetXNewtonEquationFlag()==true ) {
      for ( int j=0 ; j<this->sizeOfElement() ; j++ ) {
        if( this->element[j]->CompareNodeAddressWithFairlead( *node[i] ) ) {
          this->setJacB( n , m , j , DX , -1.0 );
        }

        if( this->element[j]->CompareNodeAddressWithAnchor( *node[i] ) ) {
          this->setJacB( n , m , j , DX , 1.0 );
        }
        n += 2;
      }
      m++;
    }


    /*
        Non-zero for (dFy)/(dH) and (dFy)/(dV)  

        Solve Y direction Newton equation
    */
    n = 0;
    if (this->node[i]->GetYNewtonEquationFlag()==true) {
      for ( int j=0 ; j<this->sizeOfElement() ; j++ ) { 
        if( this->element[j]->CompareNodeAddressWithFairlead( *node[i] ) ) {
          this->setJacB( n , m , j , DY , -1.0 ); 
        }

        if( this->element[j]->CompareNodeAddressWithAnchor( *node[i] ) ) {
          this->setJacB( n , m , j , DY , 1.0 );
        }
        n += 2;
      }
      m++; 
    }


    /*
        Non-zero for (dFz)/(dH) and (dFz)/(dV)  
        
        Solve Z direction Newton equation
    */
    n = 1;
    if (this->node[i]->GetZNewtonEquationFlag()==true) {
      for ( int j=0 ; j<this->sizeOfElement() ; j++ ) {
        if( this->element[j]->CompareNodeAddressWithFairlead( *node[i] ) ) {
          this->setJacB( n , m , j , DZ , -1.0 );
        }

        if( this->element[j]->CompareNodeAddressWithAnchor( *node[i] ) ) {
          this->setJacB( n , m , j , DZ , 1.0 );
        }
        n += 2;
      }
      m++;
    }
  }
};


// ====================================================================================================
// B_BLOCK()
// 
// The B matrix clock is based purely on the geometry of how the lines are connected to one another
//
//  J = [  A     B ]
//      [ -B^T   C ]
// 
// Get the derivatives used to populate the B block in the Jacobian.
// ====================================================================================================
double B_BLOCK::
getDeriv( )
{
  if ( partial==DX ) {
    //return 0;
    return sign*cos( this->elem->GetPsi() ); 
  } else if ( partial==DY ) {
    //return 0;
    return sign*sin( this->elem->GetPsi() ); 
  } else if ( partial==DZ ) {
    //return 0;
    return sign;
  } else {
    std::cout << "MAP : Error in B_BLOCK() function object" << std::endl; 
    return 0.0;
  };
};


void UserData::
UserJacobianEvaluations( const PetscScalar *XX )
{
  checkpoint();

  // The first step is to copy the constraint varaibles into XX (the state vector)
  for( int i=0 ; i<this->sizeOfConstraint() ; i++ ){
    this->setConstraint( i , XX[i] );
  }

  // set sum_FX, sum_FY and sum_FZ to zero
  for ( int i=0 ; i<this->sizeOfNode( ) ; i++ ){
    this->node[i]->SetSumForceToZero( );
  }
  
  // set fairlead and anchor nodes to zero. In other words, we are initializing the model. This calls
  // Node::fairlead->SetSumForceToZero();
  // Node::anchor->SetSumForceToZero();
  for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {
    this->element[i]->ResetNodes();
  }
  
  // update psi, l and h for each element being iterated
  for ( int i=0 ; i<this->sizeOfElement() ; i++ ) {
    this->element[i]->UpdateElement( *errPtr , *msgPtr );
  }

}
