/**
 * ====================================================================================================
 *                              Solving_Functions.cpp
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


#include "Element.h"  // #include "Node.h"
		      //     #include "VarType.h"
		      //         #include <boost/lexical_cast.hpp>
		      //         #include <boost/algorithm/string.hpp>
		      //         #include <string>
		      //         #include <iomanip>
		      //         #include "MAP_Message_class.h" 
		      //         #include "MAP_ErrStat_class.h"


// ====================================================================================================
// HORIZONTAL_CATENARY_EQ
//
// Horizontal catenary equation
//
// @todo : The following checks are necessary. These should be done as needed:
//         1) Specify the convergence tolerance correctly (greater than 0)
//         2) Make sure the horizontal distance between anchor and fairlead
//            is greater than 0
//         3) The fairlead cannot pass below the anchor
//         4) Unstretched length must be greater than zero
//         5) Axial stiffness must be greater than zero
//         6) The net weight of the mooring line, W, cannot be zero
//         7) check to be sure the mooring line will not double-back itself
//            and create an 'L'
// ====================================================================================================
double HORIZONTAL_CATENARY_EQ::
operator( )( )
{
  Fh = this_element->GetH( );
  Fv = this_element->GetV( );
  Lu = this_element->GetLu( );
  omega = this_element->GetOmega( );
  length = this_element->GetLength( );
  EA = this_element->GetEA( );

  //std::cout << Fh << "  " << Fv << std::endl;
  
  // check if the line is touching the seabed or not
  if ( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) ){	
    this_element->SetRestingSeabedFlag( true );
    mu = 0;
    Cb = this_element->GetCB( ); 
    Lb = Lu - Fv/omega;
	
    if ( Lb - Fh/(Cb*omega) > 0 ) { //
      mu = Lb - Fh/(Cb*omega);
    };

    // Return the equation for a chain touching the seabed
    return (Fh/omega) * boost::math::asinh( (Fv/Fh) ) 
      + ((Cb*omega)/(2*EA)) * ( -Lb*Lb + mu*(Lb - Fh/(Cb*omega)) ) + ((Fh*Lu)/(EA)) + Lb - length;
  } else {
    this_element->SetRestingSeabedFlag( false );
    
    // Free-hanging chain not in contact with the ground
    return (Fh/omega)*boost::math::asinh( Fv/Fh ) 
      - (Fh/omega)*boost::math::asinh( (Fv-omega*Lu)/Fh ) + ((Fh*Lu)/(EA)) - length;
  };// END if
};


// ====================================================================================================
// VERTICAL_CATENARY_EQ
//
// Vertical catenary equation
// 
// @todo : The following checks are necessary. These should be done as needed:
//         1) Specify the convergence tolerance correctly (greater than 0)
//         2) Make sure the horizontal distance between anchor and fairlead
//            is greater than 0
//         3) The fairlead cannot pass below the anchor
//         4) Unstretched length must be greater than zero
//         5) Axial stiffness must be greater than zero
//         6) The net weight of the mooring line, W, cannot be zero
//         7) check to be sure the mooring line will not double-back itself
//            and for an 'L'
// ====================================================================================================
double VERTICAL_CATENARY_EQ::
operator( )( )
{
  Fh = this_element->GetH( );
  Fv = this_element->GetV( );
  Lu = this_element->GetLu( );
  omega = this_element->GetOmega( );
  height = this_element->GetHeight( );
  EA = this_element->GetEA( );

  // check if the line is touching the seabed or not
  if ( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) ){	
    assert( this_element->GetRestingSeabedFlag( ) == true );
        
    // Return the equation for a chain touching the seabed
    return (Fh/omega)*(sqrt(1 + pow((Fv/Fh),2) ) - 1) + (pow(Fv,2)/(2*EA*omega)) - height;
  } else {   
    assert( this_element->GetRestingSeabedFlag( ) == false );
    
    // Return the equation for a freeling hanging chain
    return (Fh/omega)*sqrt(1 + pow( (Fv/Fh) , 2) )  
      - (Fh/omega)*sqrt( 1 + pow( ((Fv-omega*Lu)/Fh) , 2 ) ) + 1/( EA )*(Fv*Lu - (omega*Lu*Lu)/2 ) - height;
  };//END if
};


// ====================================================================================================
// SUM_F_IN_X
//
// Newton's static equilibrium sum force equation for the X direction
// ====================================================================================================
double SUM_F_IN_X::
operator( )( )
{  
  // The applied force on each each node must equation the contribution from each element
  // fairlead/anchor force applied to it. FX.value is the externally applied load
  return this_node->GetSumFX( ) - this_node->GetFX( );
};


// ====================================================================================================
// SUM_F_IN_Y
//
// Newton's static equilibrium sum force equation for the Y direction
// ====================================================================================================
double SUM_F_IN_Y::
operator( )( )
{
  // The applied force on each each node must equation the contribution from each element 
  //fairlead/anchor force  applied to it. FY.value is the externally applied load
  return this_node->GetSumFY( ) - this_node->GetFY( );
};



// ====================================================================================================
// SUM_F_IN_Z
//
// Newton's static equilibrium sum force equation for the Z direction
// ====================================================================================================
double SUM_F_IN_Z::
operator( )( )
{
  // The applied force on each each node must equation the contribution from each element fairlead/anchor 
  // force applied to it. FZ.value is the externally applied load.  We must also account for the 
  // contribution from the node mass and displaced volume/buoyancy.
  return (this_node->GetSumFZ( )) - (this_node->GetFZ( )) 
    + (this_node->GetM( )*this_node->GetGrav( )) 
    - (this_node->GetB( )*this_node->GetGrav( )*this_node->GetSeaDensity( ));
};


// ====================================================================================================
// DXDH
//
// dX/dH - partial of horizontal catenary equation with respect to the horizontal force Fh
// 
//   VFOvrHF            = VF/HF                                 Fv/Fh
//   VFOvrHF2           = VFOvrHF    //VFOvrHF                  pow( Fv/Fh , 2 ) 
//   SQRT1VFOvrHF2      = SQRT( 1.0_DbKi + VFOvrHF2      )      sqrt( 1.0 + pow( Fv/Fh , 2) )
//   VFMinWLOvrHF       = ( VF - WL )/HF;                       (Fv-omega*Lu)/Fh
//   VFMinWLOvrHF2      = VFMinWLOvrHF*VFMinWLOvrHF;            pow( (Fv-omega*Lu)/Fh , 2)
//   SQRT1VFMinWLOvrHF2 = SQRT( 1.0_DbKi + VFMinWLOvrHF2 );     sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) )
//   LOvrEA             = L  /EA                                (Lu/(EA))
//   W                  = W                                     omega
//   CBOvrEA            = CB /EA                                (Cb/(EA))
//   LMinVFOvrW         = L  - VF/W                             (Lu - Fv/omega)
//   HFOvrWEA           = HF/WEA                                (Fh/(omega*EA))
//   VFOvrWEA           = VF/WEA                                (Fv/(omega*EA))
//   HFOvrW             = HF/W                                  (Fh/omega)
//   VFMinWL            = VF - WL                               (Fv - omega*Lu)
// ====================================================================================================
double DXDH::
operator( )( )
{
  Fh     = this_element->GetH     ( );
  Fv     = this_element->GetV     ( );
  Lu     = this_element->GetLu    ( );
  omega  = this_element->GetOmega ( );
  EA     = this_element->GetEA    ( );
  Cb     = this_element->GetCB    ( );
  Lb     = Lu - Fv/omega;

//    return ( log( Fv/Fh + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) - log( (Fv-omega*Lu)/Fh + sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) ) )/omega 
//      - ( ( Fv/Fh + pow( Fv/Fh , 2 )/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( Fv/Fh + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) 
//          - ( (Fv-omega*Lu)/Fh + pow( (Fv-omega*Lu)/Fh , 2)/sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )
//          /( (Fv-omega*Lu)/Fh + sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) ) )/omega 
//      + (Lu/(EA));

  // check if the line is touching the seabed or not
  if ( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) && -Cb*(Fv-omega*Lu)<Fh ){ // if touching the seabed     
    if( Lb - Fh/(Cb*omega)>=0.0 ) { // if touching the seabed and Va > 0          
      return boost::math::asinh(Fv/Fh)/omega - ( ( (Fv/Fh) + pow(Fv/Fh , 2)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( (Fv/Fh) + sqrt( 1.0 + pow( Fv/Fh , 2) )) )/omega 
        + (Lu/(EA));
    } else { // if touching the seabed and Va = 0    
      return boost::math::asinh(Fv/Fh)/omega - ( ( (Fv/Fh) + pow(Fv/Fh , 2)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( (Fv/Fh) + sqrt( 1.0 + pow( Fv/Fh , 2) )) )/omega 
        + (Lu/(EA)) - ( (Lu - Fv/omega) - (Fh/omega)/Cb )/EA;
    }; 
  } else { // freely hanging and not touching the seabed
    return ( boost::math::asinh(Fv/Fh) - boost::math::asinh((Fv-omega*Lu)/Fh) )/omega 
      - ( ( Fv/Fh + pow( Fv/Fh , 2 )/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( Fv/Fh + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) 
          - ( (Fv-omega*Lu)/Fh + pow( (Fv-omega*Lu)/Fh , 2)/sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )
          /( (Fv-omega*Lu)/Fh + sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) ) )/omega + (Lu/(EA));    
  };


//  // check if the line is touching the seabed or not
//  if ( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) && -Cb*(Fv-omega*Lu)<Fh ){ // if touching the seabed and Va > 0    
//    return boost::math::asinh(Fv/Fh)/omega - ( ( (Fv/Fh) + pow(Fv/Fh , 2)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( (Fv/Fh) + sqrt( 1.0 + pow( Fv/Fh , 2) )) )/omega 
//      + (Lu/(EA));
////  } else if( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) && -Cb*(Fv-omega*Lu)>=Fh ) { // touching the seabed with Va == 0    
////    return boost::math::asinh(Fv/Fh)/omega - ( ( (Fv/Fh) + pow(Fv/Fh , 2)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( (Fv/Fh) + sqrt( 1.0 + pow( Fv/Fh , 2) )) )/omega 
////      + (Lu/(EA)) - ( (Lu - Fv/omega) - (Fh/omega)/Cb )/EA;    
//  } else {//if ( (Fv-omega*Lu) > 0 || omega < 0.0 )  { // freely hanging and not touching the seabed
//    return ( boost::math::asinh(Fv/Fh) - boost::math::asinh((Fv-omega*Lu)/Fh) )/omega 
//      - ( ( Fv/Fh + pow( Fv/Fh , 2 )/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( Fv/Fh + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) 
//          - ( (Fv-omega*Lu)/Fh + pow( (Fv-omega*Lu)/Fh , 2)/sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )
//          /( (Fv-omega*Lu)/Fh + sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) ) )/omega + (Lu/(EA));    
//  };
};


// ====================================================================================================
// DXDV
//
// dX/dV - partial of horizontal catenary equation with respect to the veritical force Fv
// ====================================================================================================
double DXDV::
operator( )( )
{
  Fh     = this_element->GetH     ( );
  Fv     = this_element->GetV     ( );
  Lu     = this_element->GetLu    ( );
  omega  = this_element->GetOmega ( );
  EA     = this_element->GetEA    ( );      
  Cb     = this_element->GetCB    ( );
  Lb     = Lu - Fv/omega;

//  return ( ( 1.0 + Fv/Fh /sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( Fv/Fh + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) 
//           - ( 1.0 + (Fv-omega*Lu)/Fh /sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )
//           /( (Fv-omega*Lu)/Fh + sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) ) )/omega;


  if ( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) && -Cb*(Fv-omega*Lu)<Fh ){ // if touching the seabed         
    if( Lb - Fh/(Cb*omega)>=0.0 ) { // if touching the seabed and Va > 0          
      return ( ( 1.0 + (Fv/Fh)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( (Fv/Fh) + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) )/omega + (Cb/(EA))*(Lu - Fv/omega) - 1.0/omega;    
    } else { // if touching the seabed and Va = 0    
      return ( ( 1.0 + (Fv/Fh)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( (Fv/Fh) + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) )/ omega  + (Fh/(omega*EA)) - 1.0/omega;    
    }; 
  } else { // freely hanging and not touching the seabed
    return ( ( 1.0 + Fv/Fh /sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( Fv/Fh + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) 
             - ( 1.0 + (Fv-omega*Lu)/Fh /sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )
             /( (Fv-omega*Lu)/Fh + sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) ) )/omega;
  };


//  // check if the line is touching the seabed or not
//  if ( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) && -Cb*(Fv-omega*Lu)<Fh ){ // if touching the seabed and Va > 0
//    return ( ( 1.0 + (Fv/Fh)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( (Fv/Fh) + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) )/omega + (Cb/(EA))*(Lu - Fv/omega) - 1.0/omega;    
////  } else if( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) && -Cb*(Fv-omega*Lu)>=Fh ) { // touching the seabed with Va == 0    
////    std::cout << " 2 " << std::endl;
////    return ( ( 1.0 + (Fv/Fh)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( (Fv/Fh) + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) )/ omega  + (Fh/(omega*EA)) - 1.0/omega;    
//  } else {//if ( (Fv-omega*Lu) > 0 || omega < 0.0  ) { // freely hanging and not touching the seabed
//    return ( ( 1.0 + Fv/Fh /sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( Fv/Fh + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) 
//             - ( 1.0 + (Fv-omega*Lu)/Fh /sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )
//             /( (Fv-omega*Lu)/Fh + sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) ) )/omega;
//  };
};


// ====================================================================================================
// DZDH
//
// dZ/dH - partial of vertical catenary equation with respect to the horizontal force Fh
// ====================================================================================================
double DZDH::
operator( )( )
{
  Fh     = this_element->GetH     ( );
  Fv     = this_element->GetV     ( );
  Lu     = this_element->GetLu    ( );
  omega  = this_element->GetOmega ( );
  EA     = this_element->GetEA    ( );
  Cb     = this_element->GetCB    ( );

//    return ( sqrt( 1.0 + pow( Fv/Fh , 2) ) - sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/ omega  
//      - ( pow( Fv/Fh , 2 )/sqrt( 1.0 + pow( Fv/Fh , 2) ) 
//          - pow( (Fv-omega*Lu)/Fh , 2)/sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/omega;

  if ( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) && -Cb*(Fv-omega*Lu)<Fh ){ // if touching the seabed     
    return ( sqrt( 1.0 + pow( Fv/Fh , 2) ) - 1.0 - pow(Fv/Fh , 2)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/omega;    
  } else { // freely hanging and not touching the seabed
    return ( sqrt( 1.0 + pow( Fv/Fh , 2) ) - sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/ omega  
      - ( pow( Fv/Fh , 2 )/sqrt( 1.0 + pow( Fv/Fh , 2) ) 
          - pow( (Fv-omega*Lu)/Fh , 2)/sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/omega;        
  };


//  // check if the line is touching the seabed or not
//  if ( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) && -Cb*(Fv-omega*Lu)<Fh ){ // if touching the seabed and Va > 0
//    return ( sqrt( 1.0 + pow( Fv/Fh , 2) ) - 1.0 - pow(Fv/Fh , 2)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/omega;    
////  } else if( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) && -Cb*(Fv-omega*Lu)>=Fh ) { // touching the seabed with Va == 0    
////    std::cout << " 2 " << std::endl;
////    return ( sqrt( 1.0 + pow( Fv/Fh , 2) ) - 1.0 - pow(Fv/Fh , 2)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/omega;
//  } else {//if ( (Fv-omega*Lu) > 0 || omega < 0.0 ) { // freely hanging and not touching the seabed
//    return ( sqrt( 1.0 + pow( Fv/Fh , 2) ) - sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/ omega  
//      - ( pow( Fv/Fh , 2 )/sqrt( 1.0 + pow( Fv/Fh , 2) ) 
//          - pow( (Fv-omega*Lu)/Fh , 2)/sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/omega;    
//  };
};


// ====================================================================================================
// DZDV
//
// dZ/dV - partial of vertical catenary equation with respect to the horizontal force Fv
// ====================================================================================================
double DZDV::
operator( )( )
{
  Fh     = this_element->GetH     ( );
  Fv     = this_element->GetV     ( );
  Lu     = this_element->GetLu    ( );
  omega  = this_element->GetOmega ( );
  EA     = this_element->GetEA    ( );
  Cb     = this_element->GetCB    ( );
  
//  return ( Fv/Fh/sqrt( 1.0 + pow( Fv/Fh , 2) ) - (Fv-omega*Lu)/Fh /sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/omega + (Lu/(EA));


  if ( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) && -Cb*(Fv-omega*Lu)<Fh ){ // if touching the seabed     
      return ((Fv/Fh)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/omega + (Fv/(omega*EA));            
  } else { // freely hanging and not touching the seabed
    return ( Fv/Fh/sqrt( 1.0 + pow( Fv/Fh , 2) ) - (Fv-omega*Lu)/Fh /sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/omega + (Lu/(EA));
  };


//  // check if the line is touching the seabed or not
//  if ( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) && -Cb*(Fv-omega*Lu)<Fh ){ // if touching the seabed and Va > 0
//    return ((Fv/Fh)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/omega + (Fv/(omega*EA));            
////  } else if( this_element->GetOmitContactFlag() == false && (Fv<omega*Lu) && -Cb*(Fv-omega*Lu)>=Fh ) { // touching the seabed with Va == 0    
////    std::cout << " 2 " << std::endl;
////    return ((Fv/Fh)/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/omega + (Fv/(omega*EA));
//  } else {//if ( (Fv-omega*Lu) > 0 || omega < 0.0 ) { // freely hanging and not touching the seabed
//    return ( Fv/Fh/sqrt( 1.0 + pow( Fv/Fh , 2) ) - (Fv-omega*Lu)/Fh /sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/omega + (Lu/(EA));
//  };
};
