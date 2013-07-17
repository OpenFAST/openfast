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


#include "Element.h" /**
		      *
		      * #include "Node.h"
		      *     #include "VarType.h"
		      *         #include <boost/lexical_cast.hpp>
		      *         #include <boost/algorithm/string.hpp>
		      *         #include <string>
		      *         #include <iomanip>
		      *         #include "MAP_Message.h" 
		      *         #include "MAP_ErrStat.h"
		      */

/**
 * ====================================================================================================
 * HORIZONTAL_CATENARY_EQ
 *
 * Horizontal catenary equation
 *
 * @todo : The following checks are necessary. These should be done as needed:
 *         1) Specify the convergence tolerance correctly (greater than 0)
 *         2) Make sure the horizontal distance between anchor and fairlead
 *            is greater than 0
 *         3) The fairlead cannot pass below the anchor
 *         4) Unstretched length must be greater than zero
 *         5) Axial stiffness must be greater than zero
 *         6) The net weight of the mooring line, W, cannot be zero
 *         7) check to be sure the mooring line will not double-back itself
 *            and for an 'L'
 * ====================================================================================================
 */
double HORIZONTAL_CATENARY_EQ::operator()( ){
    Fh     = this_element->getH     ( );
    Fv     = this_element->getV     ( );
    Lu     = this_element->getLu    ( );
    omega  = this_element->getOmega ( );
    length = this_element->getLength( );
    EA     = this_element->getEA   ( );

    //std::cout << Fh << "  " << Fv << std::endl;

    // check if the line is touching the seabed or not
    if ( this_element->getOMIT_CONTACT_flag() == false && (Fv<omega*Lu) ){	
        this_element->setRestingSeabed( true );
        mu = 0;
        Cb = this_element->getCB( ); 
	Lb = Lu - Fv/omega;
        
//	if( Lb > 0  && this_element->getOMIT_CONTACT_flag() == false ){
//            std::cout << "The OMIT_CONTACT flag is raised, but the cable is positively buoyant. The cable will not contact the seabed. Is the MAP input file correct? Element " << std::endl;
//        };
	
	if ( Lb - Fh/(Cb*omega) > 0 ) {
            mu = Lb - Fh/(Cb*omega);
        };// END if

	// Return the equation for a chain touching the seabed
        return (Fh/omega) * boost::math::asinh( (Fv/Fh) ) + ((Cb*omega)/(2*EA)) * ( -Lb*Lb + mu*(Lb - Fh/(Cb*omega)) ) + ((Fh*Lu)/(EA)) + Lb - length;
    }
    else {
        this_element->setRestingSeabed( false );

	// Free-hanging chain not in contact with the ground
	return (Fh/omega)*boost::math::asinh( Fv/Fh ) - (Fh/omega)*boost::math::asinh( (Fv-omega*Lu)/Fh ) + ((Fh*Lu)/(EA)) - length;
    };// END if
};


/**
 * ====================================================================================================
 * VERTICAL_CATENARY_EQ
 *
 * Vertical catenary equation
 * 
 * @todo : The following checks are necessary. These should be done as needed:
 *         1) Specify the convergence tolerance correctly (greater than 0)
 *         2) Make sure the horizontal distance between anchor and fairlead
 *            is greater than 0
 *         3) The fairlead cannot pass below the anchor
 *         4) Unstretched length must be greater than zero
 *         5) Axial stiffness must be greater than zero
 *         6) The net weight of the mooring line, W, cannot be zero
 *         7) check to be sure the mooring line will not double-back itself
 *            and for an 'L'
 * ====================================================================================================
 */
double VERTICAL_CATENARY_EQ::operator()(){
    Fh     = this_element->getH     ( );
    Fv     = this_element->getV     ( );
    Lu     = this_element->getLu    ( );
    omega  = this_element->getOmega ( );
    height = this_element->getHeight( );
    EA     = this_element->getEA   ( );

    // check if the line is touching the seabed or not
    if ( this_element->getOMIT_CONTACT_flag() == false && (Fv<omega*Lu) ){	
        assert( this_element->getRestingSeabed( ) == true );
        
        // Return the equation for a chain touching the seabed
        return (Fh/omega)*(sqrt(1 + pow((Fv/Fh),2) ) - 1) + (pow(Fv,2)/(2*EA*omega)) - height;
    }
    else{   
        assert( this_element->getRestingSeabed( ) == false );

        // Return the equation for a freeling hanging chain
        return (Fh/omega)*sqrt(1 + pow( (Fv/Fh) , 2) ) - (Fh/omega)*sqrt( 1 + pow( ((Fv-omega*Lu)/Fh) , 2 ) ) + 1/( EA )*(Fv*Lu - (omega*Lu*Lu)/2 ) - height;
    };//END if
};


/**
 * ====================================================================================================
1 * SUM_F_IN_X
 *
 * Newton's static equilibrium sum force equation for the X direction
 * ====================================================================================================
 */
double SUM_F_IN_X::operator()(){  
    // The applied force on each each node must equation the contribution from each element
    // fairlead/anchor force applied to it. FX.value is the externally applied load
    return this_node->getSum_FX( ) - this_node->getFX( );
};


/**
 * ====================================================================================================
 * SUM_F_IN_Y
 *
 * Newton's static equilibrium sum force equation for the Y direction
 * ====================================================================================================
 */
double SUM_F_IN_Y::operator()(){
    // The applied force on each each node must equation the
    // contribution from each element fairlead/anchor force 
    // applied to it. FY.value is the externally applied load
    return this_node->getSum_FY( ) - this_node->getFY( );
};


/**
 * ====================================================================================================
 * SUM_F_IN_Z
 *
 * Newton's static equilibrium sum force equation for the Z direction
 * ====================================================================================================
 */
double SUM_F_IN_Z::operator()(){
    // The applied force on each each node must equation 
    // the contribution from each element fairlead/anchor 
    // force applied to it. FZ.value is the externally applied 
    // load.  We must also account for the contribution from 
    // the node mass and displaced volume/buoyancy.
    return (this_node->getSum_FZ( )) - (this_node->getFZ( )) 
        + (this_node->getM( )*this_node->getGrav( )) 
        - (this_node->getB( )*this_node->getGrav( )*this_node->getSeaDensity( ));
};


/**
 * ====================================================================================================
 * DXDH
 *
 * dX/dH - partial of horizontal catenary equation with respect to the horizontal force Fh
 * 
 *   VFOvrHF            = VF/HF                                 Fv/Fh
 *   VFOvrHF2           = VFOvrHF     *VFOvrHF                  pow( Fv/Fh , 2 ) 
 *   SQRT1VFOvrHF2      = SQRT( 1.0_DbKi + VFOvrHF2      )      sqrt( 1.0 + pow( Fv/Fh , 2) )
 *   VFMinWLOvrHF       = ( VF - WL )/HF;                       (Fv-omega*Lu)/Fh
 *   VFMinWLOvrHF2      = VFMinWLOvrHF*VFMinWLOvrHF;            pow( (Fv-omega*Lu)/Fh , 2)
 *   SQRT1VFMinWLOvrHF2 = SQRT( 1.0_DbKi + VFMinWLOvrHF2 );     sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) )
 *   LOvrEA             = L  /EA                                (Lu/(EA))
 *   W                  = W                                     omega
 * ====================================================================================================
 */
double DXDH::operator()( ){
    Fh     = this_element->getH     ( );
    Fv     = this_element->getV     ( );
    Lu     = this_element->getLu    ( );
    omega  = this_element->getOmega ( );
    EA     = this_element->getEA   ( );
    return ( log( Fv/Fh + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) - log( (Fv-omega*Lu)/Fh + sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) ) )/omega 
        - ( ( Fv/Fh + pow( Fv/Fh , 2 )/sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( Fv/Fh + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) 
            - ( (Fv-omega*Lu)/Fh + pow( (Fv-omega*Lu)/Fh , 2)/sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/( (Fv-omega*Lu)/Fh + sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) ) )/omega 
        + (Lu/(EA));
};


/**
 * ====================================================================================================
 * DXDV
 *
 * dX/dV - partial of horizontal catenary equation with respect to the veritical force Fv
 * ====================================================================================================
 */
double DXDV::operator()( ){
    Fh     = this_element->getH     ( );
    Fv     = this_element->getV     ( );
    Lu     = this_element->getLu    ( );
    omega  = this_element->getOmega ( );
    EA     = this_element->getEA   ( );      

    return ( ( 1.0 + Fv/Fh /sqrt( 1.0 + pow( Fv/Fh , 2) ) )/( Fv/Fh + sqrt( 1.0 + pow( Fv/Fh , 2) ) ) 
             - ( 1.0 + (Fv-omega*Lu)/Fh /sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/( (Fv-omega*Lu)/Fh + sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) ) )/omega;
};


/**
 * ====================================================================================================
 * DZDH
 *
 * dZ/dH - partial of vertical catenary equation with respect to the horizontal force Fh
 * ====================================================================================================
 */
double DZDH::operator()( ){
    Fh     = this_element->getH     ( );
    Fv     = this_element->getV     ( );
    Lu     = this_element->getLu    ( );
    omega  = this_element->getOmega ( );
    EA     = this_element->getEA   ( );

    return ( sqrt( 1.0 + pow( Fv/Fh , 2) ) - sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/ omega  
      - ( pow( Fv/Fh , 2 )/sqrt( 1.0 + pow( Fv/Fh , 2) ) -  pow( (Fv-omega*Lu)/Fh , 2)/sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/omega;
};


/**
 * ====================================================================================================
 * DZDV
 *
 * dZ/dV - partial of vertical catenary equation with respect to the horizontal force Fv
 * ====================================================================================================
 */
double DZDV::operator()( ){
    Fh     = this_element->getH     ( );
    Fv     = this_element->getV     ( );
    Lu     = this_element->getLu    ( );
    omega  = this_element->getOmega ( );
    EA     = this_element->getEA   ( );

    return ( Fv/Fh/sqrt( 1.0 + pow( Fv/Fh , 2) ) - (Fv-omega*Lu)/Fh /sqrt( 1.0 + pow( (Fv-omega*Lu)/Fh , 2) ) )/omega + (Lu/(EA));
};
