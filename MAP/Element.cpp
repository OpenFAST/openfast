/**
 * ====================================================================================================
 *                              Element.cpp
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


#include "Element.h"


/**
 * ============================================================================
 * setLineProperty
 *
 * Sets a reference to a cable property in the CableLibrary. 'line_property'
 * is a pointer
 *
 * @input : element_type -- 
 * @input : T            --
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void Element::setLineProperty( const std::string              &element_type , 
			       std::vector <CableLibrary_ptr> &T            ,
			       MAP_ErrStat                    &Error        ,
			       MAP_Message                    &Msg ) {
    for( unsigned int i=0 ; i<T.size() ; i++) {

        // set the line properties. line_property is a pointer to the 
	// CableLibrary
	if ( boost::iequals( element_type , T[i]->label ) ){
	    line_property = T[i].get();
	}
    }//END for

    if (line_property == NULL ){
	// set a default value so the program does not crash
	line_property = T[0].get();
	
	// write error message to the console/glue code
	Error.set_error_key( MAP_ERROR );  
	Msg.RecordErrorToErrorList("Cannot read 'LineType' for an element. Please check for consistency between the LINE DICTIONARY and Element sections of the MAP input file."); 
    };//END if
};


/**
 * ============================================================================
 * setAnchor >
 *
 * Sets a reference to the lower node. 'anchor' is a pointer.
 *
 * @input : lower
 * @input : T
 * ============================================================================
 */
void Element::setAnchor(const int lower, const std::vector <Node_ptr> &T){
    anchor = T[lower-1].get();

    //std::cout << "Anchor node number : " << lower << std::endl;
};


/**
 * ============================================================================
 * setFairlead >
 *
 * Sets a reference to the upper node. 'fairlead' is a pointer.
 *
 * @input : T            --
 * @input : upper        -- 
 * ============================================================================
 */
void Element::setFairlead(const int upper, const std::vector <Node_ptr> &T){
    fairlead = T[upper-1].get();

    //std::cout << "Fairlead node number : " << upper << std::endl;
};


/**
 * ============================================================================
 * setH_and_V_flags >
 *
 * At initialization:
 *   -- H.is_fixed = true
 *   -- V.is_fixed = true
 *
 * H and V are iterated if:
 *   -- the fairlead or anchor connection is a 'Connect' node
 * ============================================================================
 */
void Element::setH_and_V_flags(){
    
    /**
     * =======  Set Fairlead Node flag  ======     <-------------------------------------------------------------+
     */                                                                                              //          |
                                                                                                     //          |
    // procedure for 'Connect' nodes                                                                 //          |
    if( this->fairlead->getNodeType() == Connect ){                                                  //          |
        this->setHFlagTo( false ); // have to iterate H                                              //          |
        this->setVFlagTo( false ); // have to iterate V                                              //          |
                                                                                                     //          |
        // Ensure we are not iterated FX, FY or FZ                                                   //          |
        assert( this->fairlead->FX.is_fixed == true );                                               //          |
        assert( this->fairlead->FY.is_fixed == true );                                               //          |
        assert( this->fairlead->FZ.is_fixed == true );                                               //          |
    }                                                                                                //          |
                                                                                                     //          |
    // procedure for 'Vessel' and 'Fix' nodes                                                        //          |
    else if( this->fairlead->getNodeType() == Vessel ){                                              //          |
        // For Vessel fairleads: if FX or FY value does not begin with a '#',                        //          |
        // then we are not                                                                           //          |
        if ( this->fairlead->FX.is_fixed == true || this->fairlead->FY.is_fixed == true ){           //          |
            this->setHFlagTo( true );                                                                //          |
        }                                                                                            //          |
        else{                                                                                        //          |
            this->setHFlagTo( false );                                                               //          |
        };//END if                                                                                   //          |
                                                                                                     //          |
        if ( this->fairlead->FZ.is_fixed == true ){                                                  //          |
            this->setVFlagTo( true );                                                                //          |
        }                                                                                            //          |
        else{                                                                                        //          |
            this->setVFlagTo( false );                                                               //          |
        };//END if                                                                                   //          |
    };//END if                                                                                       //   -------+
    //============== <END> =======================================================================================


    /**
     * =======  Set Anchor Node flag  ======     <---------------------------------------------------------------+
     */                                                                                              //          |
                                                                                                     //          |
     if ( this->anchor->getNodeType() == Connect ) {                                                 //          |
         // aborts the program if we enter this if                                                   //          |
         // statement and the node is attached to a                                                  //          |
         // vessel                                                                                   //          |
         assert( this->anchor->getNodeType() != Vessel );                                            //          |
                                                                                                     //          |
         // set H flag                                                                               //          |
         // if both FX and FY are fixed value, H is also fixed                                       //          |
         if ( this->fairlead->FX.is_fixed != true || this->fairlead->FY.is_fixed != true ) {         //          |
             this->setHFlagTo( false );                                                              //          |
         }                                                                                           //          |
         else {                                                                                      //          |
             this->setHFlagTo( true );                                                               //          |
         };//END if                                                                                  //          |
                                                                                                     //          |
         // set V flag                                                                               //          |
         // if FZ is fixed value, V is also fixed (assuming                                          //          |
         // M and B are fixed).                                                                      //          |
         //                                                                                          //          |
         // @todo : insert a mechanism to iterate V if M or                                          //          |
         //          B if 'is_fixed' = false                                                         //          |
         if ( this->fairlead->FZ.is_fixed != false ) {                                               //          |
             this->setVFlagTo( true );                                                               //          |
         }                                                                                           //          |
         else {                                                                                      //          |
             this->setVFlagTo( false );                                                              //          |
         };//END if                                                                                  //          |
                                                                                                     //          |
         // Ensure we are not iterating FX, FY or FZ                                                 //          |
         // becuase a connect node MUST have X, Y and                                                //          |
         // Z iterated                                                                               //          |
         assert( this->anchor->FX.is_fixed == true );                                                //          |
         assert( this->anchor->FY.is_fixed == true );                                                //          |
         assert( this->anchor->FZ.is_fixed == true );                                                //          |
     };//END if                                                                                      //   -------+
    //============== <END> =======================================================================================
    
    // @todo : is this correct? This is the exceptional case where both the
    //         fairlead and anchor nodes are defined as connect nodes
    if( this->fairlead->getNodeType() == Connect && this->fairlead->getNodeType() == Connect ){
        // have to iterate H
        this->setHFlagTo( false ); 

        // have to iterate V
        this->setVFlagTo( false ); 

        // Ensure we are not iterating FX, FY or FZ
        assert( this->fairlead->FX.is_fixed == true );
        assert( this->fairlead->FY.is_fixed == true );
        assert( this->fairlead->FZ.is_fixed == true );
    };//END if
    
    /**
     * ===============  Note  ==========================================
     *
     * An anchor node cannot be attached be attached to a vessel; hence 
     * there is no Vessel conditional check like the fairlead node  	   
     * =================================================================
     */
};


/**
 * ============================================================================
 * setLu
 *
 * Initialize Lu (unstretched cable length) variable
 *
 * @input : in           --
 * @input : i            --
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void Element::setLu( const std::string in, const int i, MAP_ErrStat &Error , MAP_Message &Msg ){
    this->Lu.name = "Lu";   // name the VarType for out file identification purposes
    this->Lu.index = i+1;   // element number, for out file identification purposes
    
    // VarType::setGenericVarType checks the input string to see if it
    //  begins with a '#' symbol. If it does, the is_fixed bool
    //  expressing is
    //   changed accordingly
    VarType::setGenericVarType( this->Lu, in, Error, Msg);
};


/**
 * ============================================================================
 * setHX
 *  
 * This variable is used solely for the purpose of writting data to the output
 * file. Hx is essentially equal to:
 *     
 *     HX.value = Hx = (H.value)*cos( psi )
 * ============================================================================
 */
void Element::setHX( ) {
    this->HX.name  = "HX";   // name the VarType for out file identification purposes
    this->HX.value = 9999.9;
};


/**
 * ============================================================================
 * setHY
 *  
 * This variable is used solely for the purpose of writting data to the output
 * file. Hx is essentially equal to:
 *     
 *     HY.value = Hy = (H.value)*sin( psi ) 
 * ============================================================================
 */
void Element::setHY( ){
    this->HY.name  = "HY";   // name the VarType for out file identification purposes
    this->HY.value = 9999.9;
};


/**
 * ============================================================================
 * setV
 *
 * Initialize V (veritical element force) variable
 *
 * @input : i            --
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void Element::setV( const int index, MAP_ErrStat &Error, MAP_Message &Msg ) {
    std::string in = "9999.9";  // initliaze variable
    this->V.name = "V";         // name the VarType for out file identification purposes
    this->V.index = index+1;    // element number, for out file identification purposes
    
    // VarType::setGenericVarType checks the input string to see if it
    // begins with a '#' symbol. If it does, the is_fixed bool expressing is
    // changed accordingly
    VarType::setGenericVarType( this->V, in, Error, Msg);
};


/**
 * ============================================================================
 * setH 
 *
 * Initialize H (horizontal element force) variable
 * 
 * @input : i            --
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 * ============================================================================
 */
void Element::setH( const int index, MAP_ErrStat &Error, MAP_Message &Msg ){
    std::string in = "9999.9";  // initliaze variable
    this->H.name = "H";         // name the VarType for out file identification purposes
    this->H.index = index+1;    // element number, for out file identification purposes
    
    // VarType::setGenericVarType checks the input string to see if it
    // begins with a '#' symbol. If it does, the is_fixed bool expressing is
    // changed accordingly
    VarType::setGenericVarType( this->H, in, Error, Msg);
};


/**
 * ============================================================================
 * getLu 
 *
 * return value of Lu (unstretched cable length)
 * ============================================================================
 */
double Element::getLu( ) const {
    return this->Lu.value;
};


/**
 * ============================================================================
 * getV
 *
 * return value of V (vertical cable force)
 *
 * @output :
 * ============================================================================
 */
double Element::getV( ) const {
    return this->V.value;
};


/**
 * ============================================================================
 * getH
 *
 * return value of H (horizontal cable force)
 *
 * @output :
 * ============================================================================
 */
double Element::getH( ) const {
    return this->H.value;
};


/**
 * ============================================================================
 * get_name
 *
 * returns the type of cable which defines this element, i.e., nylon, steel, ect
 *
 * @output :
 * ============================================================================
 */
std::string Element::getName( ) const {
    return this->line_property->label;
};


/**
 * ============================================================================
 * getLuFlag
 *
 * @output :
 * ============================================================================
 */
bool Element::getLuFlag( ) const {
    return Lu.is_fixed;
};


/**
 * ============================================================================
 * getHFlag
 *
 * @output :
 * ============================================================================
 */
bool Element::getHFlag( ) const {
    return H.is_fixed;
};


/**
 * ============================================================================
 * getVFlag
 *
 * @output : 
 * ============================================================================
 */
bool Element::getVFlag() const {
    return V.is_fixed;
};


/**
 * ============================================================================
 * getAnchorPosition >
 *
 * returns the anchor position for X, Y Z
 *
 * @input : ptr          --
 *
 * @output : double      --
 * ============================================================================
 */
double Element::getAnchorPosition( VarType Node::* ptr ) const {
    return ((*anchor).*ptr).value;
};


/**
 * ============================================================================
 * getFairleadPosition >
 *
 * returns the fairlead position for X, Y Z
 *
 * @input : ptr          --
 *
 * @output : double      -- 
 * ============================================================================
 */
double Element::getFairleadPosition( VarType Node::* ptr ) const {
    return ((*fairlead).*ptr).value;
};


/**
 * ============================================================================
 * getFairleadFlag
 *
 * Gets the 'is_fixed' boolean flag for X, Y, Z, M, B, FX, FY and FZ. Only
 * Node classs variables can be called with this expression
 *
 * @input : ptr          --
 *
 * @output : bool        -- 
 * ============================================================================
 */
bool Element::getFairleadFlag( VarType Node::* ptr ) const {
    return ((*fairlead).*ptr).is_fixed;
};


/**
 * ============================================================================
 * getAnchorFlag
 *
 * Gets the 'is_fixed' boolean flag for X, Y, Z, M, B, FX, FY and FZ. Only
 * Node classs variables can be called with this expression
 *
 * @input : ptr          --
 *
 * @output : bool        -- 
 * ============================================================================
 */
bool Element::getAnchorFlag( VarType Node::* ptr ) const {
    return ((*anchor).*ptr).is_fixed;
};


/**
 * ============================================================================
 * initializeElement
 *
 * At element initliazation, the following must be done:
 *   -- set the angle psi (orientation of xy relative to XYZ)
 *   -- set cable area and weight per unit length
 *   -- set the vertical and horizontal cable excursion
 *   -- initialize the horizontal and vertical forces applied at both the fairlead
 *      and anchor
 *
 * @input : g            --
 * @input :rho           --
 * @input :height        --
 * ============================================================================
 */
void Element::initializeElement( const double &g, const double &rho, MAP_ErrStat &err , MAP_Message &Msg ) {
    // set the angle between the cable element and the XYZ frame
    try {
        this->setPsi();
        this->initializeCableProperty( g , rho );
     } catch ( MAP_ERROR_CODE &code ) {
        std::string str = "";
        std::ostringstream S;
        S << H.index+1;
        str += boost::lexical_cast < std::string > ( code );
        str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";

        if ( code != MAP_WARNING_1 ) {
            Msg.RecordErrorToErrorList( str );
            err.set_error_key( MAP_ERROR );
            
            str.erase(2,1);
            S.str("");S.clear();
            S << ">>>> " +str +"    :    MAP error in file ";
            S << __FILE__;
            S << " on line ";
            S << __LINE__;
            S << " occurred.";
            Msg.WriteErrorToOutputFile( S.str() );
        }
    };// END try
    
    // This is the cable excursion in the element frame x and y directions
    this->set_l_and_h();
    
    // set H and V (as well as the anchor/fairlead force) 
    this->initializeHAndV( g , rho );
};


/**
 * ============================================================================
 * cleanNodes
 *
 * set the variables 'sum_FX', 'sum_FY' and 'sum_FZ' for the 'fairlead' and 
 * 'anchor' nodes to zero. This method calls Node::setSumForceToZero()
 * ============================================================================
 */
void Element::cleanNodes(){
    this->fairlead->setSumForceToZero( );    // for fairlead
    this->anchor->setSumForceToZero  ( );    // for anchor
};


/**
 * ============================================================================
 * updateElement
 *
 * The element is updated at each step in the solve process. We update:
 *   -- psi
 *   -- the cable space (l) and height (h)
 *   -- also, as psi, l and h are updated, the fairlead/anchor forces in the 
 *      element frame changes. This must be updated as well. 
 * ============================================================================
 */
void Element::updateElement( MAP_ErrStat &err , MAP_Message &Msg ) { 
    try {
        this->setPsi();
    } catch ( MAP_ERROR_CODE &code ) {
        std::string str = "";
        std::ostringstream S;
        S << H.index+1;
        str += boost::lexical_cast < std::string > ( code );
        str += "] : " + MAP_ERROR_CODE_to_string.at( code ) + S.str() + ".";

        if (code != MAP_WARNING_1 ) {
            Msg.RecordErrorToErrorList( str );
            err.set_error_key( MAP_ERROR );

            str.erase(2,1);
            S.str("");S.clear();
            S << ">>>> " +str +"    :    MAP error in file ";
            S << __FILE__;
            S << " on line ";
            S << __LINE__;
            S << " occurred.";
            Msg.WriteErrorToOutputFile( S.str() );
        };//END if
    };// END try

    this->set_l_and_h( );

    // if FX is fixed and FY is an iterated parameter:
    if ( this->fairlead->FX.is_fixed==true && this->fairlead->FY.is_fixed==false ){
        this->Hx = this->fairlead->FX.value;             
        this->Hy = tan( psi ) * this->fairlead->FX.value;
        this->H.value = sqrt( pow(this->Hx,2) + pow(this->Hy,2) );
    }

    // if FX is an iterated parameter and FY is fixed:
    else if ( this->fairlead->FX.is_fixed==false && this->fairlead->FY.is_fixed==true ){
        this->Hy = this->fairlead->FY.value;
        this->Hx = this->fairlead->FY.value/tan( psi );

        // make sure the demoniator in the above statement is not small (small
        // divisor error)

        // @change
        assert( this->psi >= 1e-3 );
        this->H.value = sqrt( pow(this->Hx,2) + pow(this->Hy,2) );
    }

    // if FX and FY are fixed and the fairlead is not a Connect node:
    else if ( this->fairlead->FX.is_fixed==true && this->fairlead->FY.is_fixed==true && this->fairlead->type != Connect ){
        this->Hx = this->fairlead->FX.value;
        this->Hy = this->fairlead->FY.value;
    }
 
    // the default case where both FX and FY are iterated
    else {
        this->Hx = this->H.value * cos( psi );
        this->Hy = this->H.value * sin( psi );   
    };// END if

    // since z and Z are in alignment, this is true at all times and a conditional statement is not needed
    this->Vz = this->V.value;

    // set anchor forces
    if ( this->is_resting_on_seabed == false ) {
        this->Ha_x = this->Hx;
        this->Ha_y = this->Hy;
        this->Va = Vz - omega*Lu.value;
    }
    else{ // element is resting on the seabed

        // @todo : This algorithm is really, really ugly. Fix this.
        this->Lb = this->Lu.value - this->V.value/this->omega;

        this->Ha_x = (this->H.value - Lb*this->getCB()*this->omega ) > 0 ? 
            (this->H.value - Lb*this->getCB()*this->omega )*cos( psi ) : 0;
        this->Ha_y = (this->H.value - Lb*this->getCB()*this->omega ) > 0 ? 
            (this->H.value - Lb*this->getCB()*this->omega )*sin( psi ) : 0;
        this->Va = 0.0;
    };// END if

    // used solely to print to output file:
    this->HX.value = this->Hx;
    this->HY.value = this->Hy;

    // add the fairlead forces to the node
    this->addForceToFairleadNode();
    this->addForceToAnchorNode();
};


/**
 * ============================================================================
 * setPsi
 * 
 * 'psi' is the angle between the element x-axis (local frame) and X-axis 
 * (global frame). This essentially produces this rotation matrix:
 *
 *           | cos(psi)  -sin(psi)  0 |
 * R(psi) =  | sin(psi)   cos(psi)  0 |
 *           |    0         0       1 |
 *
 *      1) first find psi - the angle of rotation between the element frame and the
 *      2) global reference frame
 *      3) r_j = fairlead displacement
 *      4) r_i = anchor displacement
 *      5) acos( dot( (r_j-r_i) , (u_i) ) / ( norm(r_j-r_i) ) )
 * ============================================================================
 */
void Element::setPsi( ) {
    this->norm = sqrt( pow( (this->fairlead->X.value - this->anchor->X.value) , 2) + 
                       pow( (this->fairlead->Y.value - this->anchor->Y.value) , 2) );
   
    // make sure the demoninator is not zero (or close to it) 
    // 
    // @todo : need to include an exception if the cable is perfectly vertical 
    if ( norm <= 1e-3 ) {
        throw MAP_ERROR_55;
    };// END if
    
    // @todo : not sure if this is correct; this maye create a possible conflict
    //         later on. If issues crop up with the rotation matrix, the following four 
    //         lines of code should be suspect
    if ( (this->fairlead->Y.value - this->anchor->Y.value) >= 0){
        // this simply finds the angle psi between the local and
        // global reference frames simply by evaluation trig relationships
        this->psi = acos( (this->fairlead->X.value - this->anchor->X.value)/this->norm );
    }
    else{
        this->psi = -acos( (this->fairlead->X.value - this->anchor->X.value)/this->norm );
    };// END if
};


/**
 * ============================================================================
 * set_l_and_h
 *
 * l -- the cable span is equal to the norm of the difference between the 
 *      fairlead and anchor displacement vectors (x and y components only) 
 *
 * h -- the z and Z axis are in alignment; so the cable height is equal to the 
 *      difference between the anchor and fairlead (z component only)
 * ============================================================================
 */
inline void Element::set_l_and_h( ) {
    // set l (horizontal displacement between fairlead and anchor) and h (vertical
    // displacement between fairlead and anchor)
    l  = sqrt( pow( (this->fairlead->X.value - this->anchor->X.value) , 2) + 
               pow( (this->fairlead->Y.value - this->anchor->Y.value) , 2) );
    
    h = this->fairlead->Z.value - this->anchor->Z.value;
};


/**
 * ============================================================================
 * initializeCableProperty
 *
 * @input : g            --
 * @input : rho          -- 
 * ============================================================================
 */
void Element::initializeCableProperty( const double &g , const double &rho ) {
    A = MAP_CONST::PI*pow( (this->line_property->Diam.value/2) , 2 );  // cross-section area
    omega = g*(this->line_property->MassDenInAir.value) - g*A*rho;       // weight per length

    // make sure the cable is not neutrally buoyant
    // if the cable sensity matches the sea density, then the program will 
    // terminate
    if ( fabs( omega ) <= 1 ){    
        throw MAP_ERROR_56;
    };//END if

    // if OMIT_CONTACT and cable density is less than seawater, throw a warning
    if( omega > 0 && OMIT_CONTACT_flag == true ) {
        throw MAP_WARNING_1;
    };// END if
};


/**
 * ============================================================================
 * initializeHAndV
 *
 * If H and V are not user supplied initial guesses or fixed values, then the 
 * process to estimate H and V is invoked using the method given in:
 *
 * @article{peyrot1979,
 *     title={{Analysis of Cable Structuresx}},
 *     author={Peyrot, A. H., and Goulois, A. M.},
 *     journal={Computers and Structures},
 *     volume={10},
 *     pages={805--813},
 *     year={1979},
 *     publisher={Pergamon Press, Ltd.}
 * }
 *
 * @input : g            --
 * @input : rho          --
 * ============================================================================
 */
void Element::initializeHAndV( const double &g , const double &rho) {

    // initialize lambda
    double lambda = 9999.9;

    if ( l==0 ){
        lambda = 1000000;
    }
    else if ( sqrt( pow(l,2) + pow(h,2) ) >= this->Lu.value ){
        lambda = 0.2;
    }
    else{
        lambda = sqrt(3*( (pow(this->Lu.value,2) - pow(h,2))/pow(l,2) - 1));
    };//END if

    assert( lambda != 9999.9 );

    // set up temporary variables to store the estimate forces from Peyrot
    // & Goulois
    const double tempH = abs( omega*l / ( 2*lambda ) );
    const double tempV = ( omega/2 ) * ( (h/tanh( lambda ) ) + Lu.value );

    // fairlead node forces
    if (this->fairlead->type != Connect){

        /**
         * =======  Set X direction fairlead force  ======     <--------------------------+
         */                                                                   //          |
                                                                              //          |
        if ( this->fairlead->FX.is_fixed==true ){                             //          |
            Hx = this->fairlead->FX.value;                                    //          |
        }                                                                     //          |
        else{                                                                 //          |
            Hx = tempH * cos( psi );                                          //          |
        };// END if                                                           //----------+
        //============== <END> ============================================================

        /**
         * =======  Set Y direction fairlead force  ======     <--------------------------+
         */                                                                   //          |
                                                                              //          |
        if ( this->fairlead->FY.is_fixed==true ){                             //          |
            Hy = this->fairlead->FY.value;                                    //          |
        }                                                                     //          |
        else{                                                                 //          |
            Hy = tempH * sin( psi );                                          //          |
        };// END if                                                           //----------+
        //============== <END> ============================================================

        /**
         * =======  Set Z direction fairlead force  ======     <--------------------------+
         */                                                                   //          |
                                                                              //          |
        if ( this->fairlead->FZ.is_fixed==true ){ // if the FZ has a fixed value          | 
            Vz = this->fairlead->FZ.value;                                    //          |
            Vz = this->fairlead->FZ.value - this->fairlead->M.value*g         //          |
                + this->fairlead->B.value*( g*rho );                          //          |
        }                                                                     //          |
        else{                                                                 //          |
            Vz = tempV;                                                       //          |
        };// END if                                                           //----------+
        //============== <END> ============================================================
    }
    else{ // this is the default for a connect node
        // This is H and V in the local element frame
        Hx = tempH * cos( psi );
        Hy = tempH * sin( psi );
        Vz = tempV;
    };// END if

    // the value for H (in the element frame) is equal to the vector 
    // sum of Hx and Hy
    this->H.value = sqrt( pow(Hx,2) + pow(Hy,2) );
    this->V.value = Vz;

    // anchor node forces
    if (this->anchor->type != Connect){

        /**
         * =======  Set Z direction anchor force  ======     <----------------------------+
         */                                                                   //          |
                                                                              //          |
        // if the FX has a fixed value                                        //          |
        if ( this->anchor->FX.is_fixed==true ){                               //          |
            Ha_x = this->anchor->FX.value;                                    //          |
        }                                                                     //          |
        else{                                                                 //          |
            Ha_x = H.value * cos( psi );                                      //          |
        };// END if                                                           //----------+
        //============== <END> ============================================================
        
        /**
         * =======  Set Z direction anchor force  ======     <----------------------------+
         */                                                                   //          |
                                                                              //          |
        // if the FY has a fixed value                                        //          |
        if ( this->anchor->FY.is_fixed==true ) {                              //          |
            Ha_y = this->anchor->FY.value;                                    //          |
        }                                                                     //          |
        else{                                                                 //          |
            Ha_y = H.value * sin( psi );                                      //          |
        };// END if                                                           //----------+
        //============== <END> ============================================================

        /**
         * =======  Set Z direction anchor force  ======     <----------------------------+
         */                                                                   //          |
                                                                              //          |
        // if the FZ has a fixed value                                        //          |
        if ( this->anchor->FZ.is_fixed==true ) {                              //          |
            Va = this->anchor->FZ.value;                                      //          |
            Va = this->anchor->FZ.value - this->anchor->M.value*g             //          |
                + this->anchor->B.value*( g*rho );                            //          |
        }                                                                     //          |
        else{                                                                 //          |
            Va = Vz - omega*Lu.value;                                         //          |
        };// END if                                                           //----------+
        //============== <END> ============================================================
    } 
    else{ // this is the default for a connect node
        // This is H and V in the local element frame
        Ha_x = H.value * cos( psi );
        Ha_y = H.value * sin( psi );
        Va = Vz - omega*Lu.value;
    };//END if

    /**
     * ===============  Update Node Force  =============================
     *
     * The last step is to update the FX, FY and FZ. Fi is the force 
     * applied to the node to hold it in static equilibrium. This is 
     * essentially Hx, Hy and V (minus the effects of the buoyancy 
     * module and added node weight)
     *
     * ======  !!!!!!  ====== 
     *
     * THE NEXT LINES CHANGE VALUES SEEN IN MAP_OtherStateType_class
     *
     * ======  !!!!!!  ====== 
     *
     * An anchor node cannot be attached be attached to a vessel; hence 
     * there is no Vessel conditional check like the fairlead node  	   
     *
     * =================================================================
     */

    this->addForceToFairleadNode();
    this->addForceToAnchorNode();
};


/**
 * ============================================================================
 * addForceToFairleadNode
 * ============================================================================
 */
inline void Element::addForceToFairleadNode(){  
    // now sum this force in with the other applied node forces
    this->fairlead->addToSumFX( Hx );
    this->fairlead->addToSumFY( Hy );
    this->fairlead->addToSumFZ( Vz );
};


/**
 * ============================================================================
 * addForceToAnchorNode
 *
 * We add the negative of the element anchor force to the node (Node::sum_FX,
 * Node::sum_FY, Node::sum_FY)
 * ============================================================================
 */
inline void Element::addForceToAnchorNode(){
    // the anchor foces are negative because they act opposite of the 
    // force applied at the fairlead
    this->anchor->addToSumFX( -Ha_x );
    this->anchor->addToSumFY( -Ha_y );
    this->anchor->addToSumFZ( -Va );
};


/**
 * ============================================================================
 * elementPlotPoints
 *
 * Allows the cable to be plotted in Python using the matplotlib library
 *
 * @input  : X           -- X global position of element point [meters]
 * @input  : Y           -- Y global position of element point [meters]
 * @input  : Z           -- Z global position of element point [meters]
 * @input  : Msg         -- Error message status
 *
 * @output : integer     -- error code. If '1', no error. If '-1', an error is 
 *                          reported. 
 * ============================================================================
 */
bool Element::elementPlotPoints( std::string &X , 
                                 std::string &Y ,
                                 std::string &Z ,
                                 MAP_Message &Msg) {
   
    int    N     = 100;                              // number of points being plotted
    double n     = boost::lexical_cast<double>( N ); // number of points being plotted 
    double Fh    = this->H.value;                    // horizontal fairlead force
    double Fv    = this->V.value;                    // vertical fairlead force
    double w     = this->omega;			     // mass per length 
    double lu    = this->Lu.value;		     // unstretched length
    //double a     = this->A;			     // cross-sectional area
    //double e     = this->getE();  		     // young's modulus
    double ea    = this->getEA();                    // element stiffness
    double angle = this->psi;                        // angles between element and global axis
    double cb    = this->getCB();                    // cable/sea bed friction coefficient

    double x = 9999.9;                               // X global position of point on element
    double y = 9999.9;                               // Y global position of point on element
    double z = 9999.9;                               // Z global position of point on element
    double S = 9999.9;                               // position along cable element
    
    std::string temp_X = "";     
    std::string temp_Y = "";      
    std::string temp_Z = "";  

    for ( int s=0 ; s<n+1 ; s++ ){

        S = (s/n)*lu;

        // If the cable is not resting on the seabed, we use the classic catenary equation
        // for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
        // the modified version as done in the FAST wind turbine program. 
        //
        // @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
        //        Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.

        /**
         * =======  Catenary equations for a hanging chain  ======     <--------------------------------------------+
         */                                                                                             //          |
                                                                                                        //          |
        if (this->is_resting_on_seabed == false ){ // Element is NOT resting on sea floor               //          |
                                                                                                        //          |
            // X position of element in global coordinates                                              //          |
            x = ((Fh/w)*boost::math::asinh( Fv/Fh )                                                     //          |
                 - (Fh/w)*boost::math::asinh( (Fv-S*w)/Fh )                                             //          |
                 + (Fh*S)/(ea))*cos( angle ) - this->fairlead->X.value;                                //          |
                                                                                                        //          |
            // Y position of element in global coordinates                                              //          |
            y = (  (Fh/w)*boost::math::asinh( Fv/Fh )                                                   //          |
                   - (Fh/w)*boost::math::asinh( (Fv-S*w)/Fh )                                           //          |
                   + (Fh*S)/(ea))*sin( angle ) - this->fairlead->Y.value;                              //          |
                                                                                                        //          |
            // Z position of element in global coordinates                                              //          |
            z =  (Fh/w)*( sqrt( 1+pow(Fv/Fh,2) ) - sqrt( 1+pow((Fv-w*S)/Fh,2) ) )                       //          |
                + (1/(ea))*(Fv*S+w*S*S/2) - this->fairlead->Z.value;                                   //          |
        }// END if                                                                                      //----------+
        //============== <END> =======================================================================================


        /**
         * =======  Catenary equations for a cable in contact with the sea bed  ======     <-------------------------+
         */                                                                                              //          |
                                                                                                         //          |
        else{ // Element is resting on sea floor                                                         //          |
                                                                                                         //          |
            Lb     = lu - (Fv/w);                                                                        //          |
            double lambda = 0.0;                                                                         //          |
                                                                                                         //          |
            // Find the right equation for X(s) and Y(s) to use for plotting the element                 //          |
            // resting on the seabed                                                                     //          |
            if ( 0<=S && S<=(Lb-Fh/(cb*w)) ) { // for 0 <= s <= Lb - H/(Cb*w)                            //          |
                // X position of element in global coordinates                                           //          |
                x = -(S) * cos( angle ) - this->anchor->X.value;                                         //          |
                                                                                                         //          |
                // Y position of element in global coordinates                                           //          |
                y = -(S) * sin( angle ) - this->anchor->Y.value;                                         //          |
            }                                                                                            //          |
            else if ( (Lb-Fh/(cb/w))<S && S<=Lb )  { // for Lb - H/(Cb*w) < s <= Lb                      //          |
                                                                                                         //          |
                lambda = (Lb-Fh/(cb*w))>0 ? (Lb-Fh/(cb*w)) : 0;                                          //          |
                                                                                                         //          |
                // X position of element in global coordinates                                           //          |
                x = -( S + ((cb*w)/(2*ea)) * (S*S - 2*(Lb-Fh/(cb*w))*S + (Lb- Fh/(cb*w))*lambda ) )     //          |
                    * cos( angle ) - this->anchor->X.value;                                              //          |
                                                                                                         //          |
                // Y position of element in global coordinates                                           //          |
                y = -( S + ((cb*w)/(2*ea)) * (S*S - 2*(Lb-Fh/(cb*w))*S + (Lb- Fh/(cb*w))*lambda ) )     //          |
                    * sin( angle ) - this->anchor->Y.value;                                              //          |
            }                                                                                            //          |
            else { // for Lb < s <= L                                                                    //          |
                                                                                                         //          |
                lambda = (Lb-Fh/(cb*w))>0 ? (Lb-Fh/(cb*w)) : 0;                                          //          |
                                                                                                         //          |
                // X position of element in global coordinates                                           //          |
                x = -( Lb + (Fh/w)*boost::math::asinh( (w*(S-Lb))/Fh ) ) * cos( angle )                  //          |
                    + ( ((Fh*S)/(ea)) + ((cb*w)/(2*ea))*( -Lb*Lb + (Lb-Fh/(cb*w))*lambda ) )           //          |
                    * cos( angle ) - this->anchor->X.value;                                              //          |
                                                                                                         //          |
                // Y position of element in global coordinates                                           //          |
                y = -( Lb + (Fh/w)*boost::math::asinh( (w*(S-Lb))/Fh ) ) * sin( angle )                  //          |
                    + ( ((Fh*S)/(ea)) + ((cb*w)/(2*ea))*( -Lb*Lb + (Lb-Fh/(cb*w))*lambda ) )           //          |
                    * sin( angle ) - this->anchor->Y.value;                                              //          |
                                                                                                         //          |
            }// END if                                                                                   //          |
                                                                                                         //          |
            // Find the right equation for Z(s) to use for plotting the element                          //          |
            // resting on the seabed                                                                     //          |
            if ( 0<=S && S<=Lb ) {                                                                       //          |
                // Z position of element in global coordinates                                           //          |
                z = 0 - this->anchor->Z.value;                                                           //          |
            }                                                                                            //          |
            else{                                                                                        //          |
                // Z position of element in global coordinates                                           //          |
                z = -( (Fh/w)*( sqrt(1 + pow((w*(S-Lb)/Fh),2) ) - 1 ) + ((w*pow((S-Lb),2))/(2*ea)) )    //          |
                    - this->anchor->Z.value;                                                             //          |
            }//END if                                                                                    //          |
        };//END if                                                                                       //----------+
        //============== <END> =======================================================================================



        // Now check to make sure the valid numbers are calculated for the x, y and z 
        // direction. The following lines of code convert real numbers into string arguments
        // for plotting purposes. 
        //
        // If -1 is returned, an exception is raised by boost::lexical_cast
        
        // X direction
        try{	    
            X += boost::lexical_cast<std::string>( -x );
            if ( s != n ) X += " , ";
        }
        catch ( boost::bad_lexical_cast const& ) {
            Msg.RecordWarningToWarningList( "'X' displacement is not a valid real number (possible NaN)." );
            return false;
        };//END try

        // Y direction
        try{
            Y += boost::lexical_cast<std::string>( -y );
            if ( s != n ) Y += " , ";
        }
        catch ( boost::bad_lexical_cast const& ) {
            Msg.RecordWarningToWarningList( "'Y' displacement is not a valid real number (possible NaN)." );
            return false;
        };//END try

        // Z direction
        try{
            Z += boost::lexical_cast<std::string>( -z );
            if ( s != n ) Z += " , ";
        }
        catch ( boost::bad_lexical_cast const& ) {
            Msg.RecordWarningToWarningList( "'Z' displacement is not a valid real number (possible NaN)." );
            return false;
        };//END try
    };
   
    return true;
};


/**
 * ============================================================================
 * getLineTension
 *
 * Prints the tension in the fairlead to the MAP map.out file
 *
 * @output : string           -- tension in for 10 points along the cable 
 *                               [Newtons]
 *
 * @todo :                    -- right now, the code is hard coded to evaluate 
 *                               10 points along the line. In the future, we 
 *                               may want to change this so an arbitrary
 *                               number of points are evaluated
 * ============================================================================
 */
std::string Element::getLineTension( ) {

    std::string        out         = ""; // initlize the string we are writting for the outputs
    std::string        temp        = "";
    std::ostringstream S;
    double             length      = 9999.9;
    double             temp_double = 9999.9;
    
    S << std::fixed << std::setprecision(3);
    
    double Fh    = this->H.value;                    // horizontal fairlead force
    double Fv    = this->V.value;                    // vertical fairlead force
    double w     = this->omega;			     // mass per length 
    double cb    = this->getCB();                    // cable/sea bed friction coefficient
    double lu    = this->Lu.value;		     // unstretched length
    double Lb    = lu - (Fv/w);			     
    
    for (int i=0 ; i<10 ; i++){ // i<10 since we are only evaluating 10 points along cable

        // @todo : right now, the code is 
        //         hard coded to evaluate
        //         10 points along the line.
        //         In the future, we may want
        //         to change this so an arbitrary
        //         number of points are evaluated

        length = (static_cast<double>(i)/9)*(this->Lu.value); 

        /**
         * =======  Case for cable hanging in air  ======     <------------------------------------------------------+
         */                                                                                              //          |
                                                                                                         //          |
        if ( this->is_resting_on_seabed==false  ) { // Element is NOT resting on sea floor               //          |
            S << sqrt( pow((this->H.value),2) + pow((this->Va+this->omega*length),2) );                  //          |
            temp += S.str();                                                                             //          |
            S.str("");S.clear();                                                                         //          |
                                                                                                         //          |
            // make sure the string array is long enough to put at least one white space                 //          |
            // at the tail; otherwise, the do-while loop will result in a program crash                  //          |
            assert( temp.size() < 14 );                                                                  //          |
                                                                                                         //          |
            // fill string with white spaces to bring it to the correct length                           //          |
            do {                                                                                         //          |
                temp += " ";                                                                             //          |
            } while( temp.size()<15 ); // @todo : 15 is arbitrarily chosen.                              //          |
                                                                                                         //          |
            out += temp;                                                                                 //          |
            temp = "";                                                                                   //          |
        }                                                                                                //----------+
        //============== <END> =======================================================================================

        /**
         * =======  Case for cable in contact with the sea bed  ======     <-----------------------------------------+
         */                                                                                              //          |
                                                                                                         //          |
        else {                                                                                           //          |
            if ( 0<=length && length<=Lb ){ // for 0 <= s <= Lb                                          //          |
                temp_double = (Fh+cb*w*(length-Lb))>0 ? (Fh+cb*w*(length-Lb)) : 0;                       //          |
                S << temp_double;                                                                        //          |
            }                                                                                            //          |
            else {// for Lb < s <= L                                                                     //          |
                S << sqrt( Fh*Fh + pow( (w*(length-Lb)) ,2) );                                           //          |
            };//END if                                                                                   //          |
                                                                                                         //          |
            temp += S.str();                                                                             //          |
            S.str("");S.clear();                                                                         //          |
                                                                                                         //          |
            // make sure the string array is long enough to put at least one white space                 //          |
            // at the tail; otherwise, the do-while loop will result in a program crash                  //          |
            assert( temp.size() < 14 );                                                                  //          |
                                                                                                         //          |
            // fill string with white spaces to bring it to the correct length                           //          |
            do {                                                                                         //          |
                temp += " ";                                                                             //          |
            } while( temp.size()<15 ); // @todo : 15 is arbitrarily chosen.                              //          |
                                                                                                         //          |
            out += temp;                                                                                 //          |
            temp = "";                                                                                   //          |
        };//END if                                                                                       //----------+
        //============== <END> =======================================================================================
    };// END for

    return out;
};


/**
 * ============================================================================
 * getLineTensionHeader 
 *
 * Prints the heades "T[1] T[2] ... T[N]" header for the element in the MAP
 * map.out file
 *
 * @input  : int i             -- the element number
 *
 * @output : string            -- "T[1] T[2] ... T[N]"
 * ============================================================================
 */
std::string Element::getLineTensionHeader( const int i ) {
    
    std::string        out  = ""; // initialize the string we are writing for the outputs
    std::string        temp = "";
    std::ostringstream S;
    
    S << std::fixed << std::setprecision(1);

    for (int j=0 ; j<10 ; j++){

        // create the string we are writting as an output
        temp += "T[";
        S << i+1;
        temp += S.str();
        S.str("");S.clear();
        temp += "](";

        // @todo : right now, the code is 
        //         hard coded to evaluate
        //         10 points along the line.
        //         In the future, we may want
        //         to change this so an arbitrary
        //         number of points are evaluated
        S <<  (static_cast<double>(j)/9)*(this->Lu.value);
 
        temp += S.str();
        temp += ")";
        S.str("");S.clear();

        // make sure the string array is long enough to put at least one white space 
        // at the tail; otherwise, the do-while loop will result in a program crash
        assert( temp.size() < 14 );

        // fill string with white spaces to bring it to the correct length
        do {
            temp += " ";
        } while( temp.size()<15 ); // @note : 15 is arbitrarily chosen. Fix this ? 

        out += temp;
        temp = "";
    };// END for 
    
     return out;
};


/**
 * ============================================================================
 * compareNodeAddressWithFairlead
 * ============================================================================
 */
bool Element::compareNodeAddressWithFairlead( const Node &ref ) const { 
   
    if ( &ref==this->fairlead ) return true;
    else return false;
};


/**
 * ============================================================================
 * compareNodeAddressWithAnchor
 * ============================================================================
 */
bool Element::compareNodeAddressWithAnchor( const Node &ref ) const {    
    if ( &ref==this->anchor ) return true;
    else return false;
};
