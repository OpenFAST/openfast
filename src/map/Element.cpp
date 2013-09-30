/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   Element.cpp
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


#include "Element.h"


/**
 * SetLineProperty
 *
 * Sets a reference to a cable property in the CableLibrary. 'line_property'
 * is a pointer
 *
 * @input : element_type -- 
 * @input : T            --
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 */
void Element::
SetLineProperty( const std::string              &element_type , 
                 std::vector <CableLibrary_ptr> &library      ,
                 MAP_ErrStat_class              &err          ,
                 MAP_Message_class              &msg          ) 
{
  for( unsigned int i=0 ; i<library.size() ; i++) {
    
    // set the line properties. line_property is a pointer to the CableLibrary
    if ( boost::iequals( element_type , library[i]->label ) ){
      this->line_property = library[i].get();
    }
  }
  
  if ( line_property == NULL ) {
    // set a default value so the program does not crash, but so that it can end gracefully. 
    this->line_property = library[0].get();

    std::string str = "";
    MAPSetUniversalErrorStat( MAP_ERROR_77 , str, err, msg );    
  }
};


/**
 * SetAnchor
 *
 * Sets a reference to the lower node. 'anchor' is a pointer.
 *
 * @input : lower
 * @input : T
 */
void Element::
SetAnchor(const int                    lower            , 
          const std::vector <Node_ptr> &node_library    )
{
  anchor = node_library[lower-1].get(); // get the address of a node
};


/** 
 * SetFairlead 
 *
 * Sets a reference to the upper node. 'fairlead' is a pointer.
 *
 * @input : T            --
 * @input : upper        -- 
 */
void Element::
SetFairlead(const int upper                            , 
            const std::vector <Node_ptr> &node_library )
{
  fairlead = node_library[upper-1].get();
};


/**
 * SetH_and_V_flags >
 *
 * At initialization:
 *   -- H.is_fixed = true
 *   -- V.is_fixed = true
 *
 * H and V are iterated if:
 *   -- the fairlead or anchor connection is a 'Connect' node
 */
void Element::
SetH_and_V_flags( )
{
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                            N O T E

     An anchor node cannot be attached be attached to a vessel; hence 
     there is no Vessel conditional check like the fairlead node  	   
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  */

  // Set Fairlead Node flag  
  if( this->fairlead->GetNodeType() == Connect ){ 
    this->SetHFlagTo( false ); // __MUST__ have to iterate H for this element
    this->SetVFlagTo( false ); // __MUST__ have to iterate V for this element

    // Ensure we are not iterated FX, FY or FZ     
    assert( this->fairlead->FX.is_fixed == true ); 
    assert( this->fairlead->FY.is_fixed == true ); 
    assert( this->fairlead->FZ.is_fixed == true ); 
  } else if( this->fairlead->GetNodeType() == Vessel ) { 
    
    // For Vessel fairleads: if FX or FY value does not begin with a '#', then we are not iterating it
    if ( this->fairlead->FX.is_fixed == true || this->fairlead->FY.is_fixed == true ){ 
      this->SetHFlagTo( true );                                                        
    } else {                                                                              
      this->SetHFlagTo( false );                                                       
    }

    if ( this->fairlead->FZ.is_fixed == true ){                                        
      this->SetVFlagTo( true );                                                        
    } else {                                                                              
      this->SetVFlagTo( false );                                                       
    }
  }


  // Set Anchor Node flag <---------------------------------------------------------------+
  if ( this->anchor->GetNodeType() == Connect ) {     
    // set H flag if both FX and FY are fixed value, H is also fixed 
    if ( this->fairlead->FX.is_fixed != true || this->fairlead->FY.is_fixed != true ) { 
      this->SetHFlagTo( false );                                                        
    } else {                                                            
      this->SetHFlagTo( true );                                                         
    };
    
    // set V flag if FZ is fixed value, V is also fixed (assuming M and B are fixed).
    // @todo : insert a mechanism to iterate V if M or B if 'is_fixed' = false                                                 
    if ( this->fairlead->FZ.is_fixed != false ) {                                       
      this->SetVFlagTo( true );                                                         
    } else {                                                                              
      this->SetVFlagTo( false );                                                        
    }

    // Ensure we are not iterating FX, FY or FZ becuase a connect node MUST have X, Y and Z iterated
    assert( this->anchor->FX.is_fixed == true );                                        
    assert( this->anchor->FY.is_fixed == true );                                        
    assert( this->anchor->FZ.is_fixed == true );                                        
  }
    
  // @todo : is this correct? This is the exceptional case where both the
  //         fairlead and anchor nodes are defined as connect nodes
  if( this->fairlead->GetNodeType() == Connect && this->fairlead->GetNodeType() == Connect ){    
    this->SetHFlagTo( false ); // have to iterate H
    this->SetVFlagTo( false ); // have to iterate V

    // Ensure we are not iterating FX, FY or FZ
    assert( this->fairlead->FX.is_fixed == true );
    assert( this->fairlead->FY.is_fixed == true );
    assert( this->fairlead->FZ.is_fixed == true );
  }
};


/**
 * SetLu
 *
 * Initialize Lu (unstretched cable length) variable
 *
 * @input : in           --
 * @input : i            --
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 */
void Element::
SetLu( const std::string varStr , 
       const int         i      , 
       MAP_ErrStat_class &err   , 
       MAP_Message_class &msg   )
{
    this->Lu.name = "Lu";   // name the VarType for out file identification purposes
    this->Lu.index = i+1;   // element number, for out file identification purposes
    
    // VarType::SetGenericVarType checks the input string to see if it
    // begins with a '#' symbol. If it does, the is_fixed bool
    // expressing is changed accordingly
    VarType::SetGenericVarType( this->Lu, varStr, err, msg );
};


/**
 * SetHX
 *  
 * This variable is used solely for the purpose of writting data to the output
 * file. Hx is essentially equal to:
 *     
 *     HX.value = Hx = (H.value)*cos( psi )
 */
void Element::
SetHX( ) 
{
  this->HX.name  = "HX";   // name the VarType for out file identification purposes
  this->HX.value = 9999.9;
};


/**
 * SetHY
 *  
 * This variable is used solely for the purpose of writting data to the output
 * file. Hx is essentially equal to:
 *     
 *     HY.value = Hy = (H.value)*sin( psi ) 
 */
void Element::
SetHY( )
{
  this->HY.name  = "HY";   // name the VarType for out file identification purposes
  this->HY.value = 9999.9;
};


/**
 * SetV
 *
 * Initialize V (veritical element force) variable
 *
 * @input : i            --
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 */
void Element::
SetV( const int         index , 
      MAP_ErrStat_class &err  , 
      MAP_Message_class &msg  ) 
{
  std::string varStr = "9999.9";  // initliaze variable
  this->V.name = "V";         // name the VarType for out file identification purposes
  this->V.index = index+1;    // element number, for out file identification purposes
  
  // VarType::SetGenericVarType checks the input string to see if it begins with
  // a '#' symbol. If it does, the is_fixed bool expressing is changed accordingly
  VarType::SetGenericVarType( this->V, varStr, err, msg);
};


/**
 * SetH 
 *
 * Initialize H (horizontal element force) variable
 * 
 * @input : i            --
 * @input : Msg          -- Error message status
 * @input : Error        -- Error code
 */
void Element::
SetH( const int         index  , 
      MAP_ErrStat_class &err ,
      MAP_Message_class &msg )
{
  std::string varStr = "999.9";  // initliaze variable
  this->H.name = "H";            // name the VarType for out file identification purposes
  this->H.index = index+1;       // element number, for out file identification purposes
  
  // VarType::SetGenericVarType checks the input string to see if it begins with
  // a '#' symbol. If it does, the is_fixed bool expressing is changed accordingly
  VarType::SetGenericVarType( this->H, varStr, err, msg );
};


/**
 * GetLu 
 *
 * return value of Lu (unstretched cable length)
 */
double Element::
GetLu( ) const 
{
  return this->Lu.value;
};


/**
 * GetV
 *
 * return value of V (vertical cable force)
 *
 * @output :
 */
double Element::
GetV( ) const 
{
  return this->V.value;
};


/**
 * GetH
 *
 * return value of H (horizontal cable force)
 *
 * @output :
 */
double Element::
GetH( ) const 
{
  return this->H.value;
};


/**
 * get_name
 *
 * returns the type of cable which defines this element, i.e., nylon, steel, ect
 *
 * @output :
 */
std::string Element::
GetElementName( ) const 
{
  return this->line_property->label;
};


/**
 * GetLuFlag
 *
 * @output :
 */
bool Element::
GetLuFlag( ) const
{
  return Lu.is_fixed;
};


/**
 * GetHFlag
 *
 * @output :
 */
bool Element::
GetHFlag( ) const 
{
  return H.is_fixed;
};


/**
 * GetVFlag
 *
 * @output : 
 */
bool Element::
GetVFlag() const
{
  return V.is_fixed;
};


/**
 * GetAnchorPosition
 *
 * returns the anchor position for X, Y Z
 *
 * @input : ptr          --
 *
 * @output : double      --
 */
double Element::
GetAnchorPosition( VarType Node::* ptr ) const
{
  return ((*anchor).*ptr).value;
};


/**
 * GetFairleadPostion 
 *
 * returns the fairlead position for X, Y Z
 *
 * @input : ptr          --
 *
 * @output : double      -- 
 */
double Element::
GetFairleadPostion( VarType Node::* ptr ) const 
{
  return ((*fairlead).*ptr).value;
};


/**
 * GetFairleadFlag
 *
 * Gets the 'is_fixed' boolean flag for X, Y, Z, M, B, FX, FY and FZ. Only
 * Node classs variables can be called with this expression
 *
 * @input : ptr          --
 *
 * @output : bool        -- 
 */
bool Element::
GetFairleadFlag( VarType Node::* ptr ) const 
{
  return ((*fairlead).*ptr).is_fixed;
};


/**
 * GetAnchorFlag
 *
 * Gets the 'is_fixed' boolean flag for X, Y, Z, M, B, FX, FY and FZ. Only
 * Node classs variables can be called with this expression
 *
 * @input : ptr          --
 *
 * @output : bool        -- 
 */
bool Element::
GetAnchorFlag( VarType Node::* ptr ) const 
{
  return ((*anchor).*ptr).is_fixed;
};


// ============================================================================
// InitializeElement
//
// At element initliazation, the following must be done:
//   -- set the angle psi, orientation of xy relative to XYZ. (not needed?)
//   -- set cable area and weight per unit length
//   -- set the vertical and horizontal cable excursion
//   -- initialize the horizontal and vertical forces applied at both the fairlead
//      and anchor
//
// @input : g            --
// @input :rho           --
// @input :height        --
// ============================================================================
void Element::
InitializeElement( const double &g       , 
                   const double &rho_sea , 
                   MAP_ErrStat_class  &err     , 
                   MAP_Message_class  &msg     ) 
{
  try {
    InitializeCableProperty( g , rho_sea ); // throws MAP_ERROR_56 and/or MAP_WARNING_1
  } catch ( MAP_ERROR_CODE &code ) { 
    std::ostringstream S;
    std::string str = "";
    S << H.index;
    str = S.str() + ".";

    // @todo : this throw a MAP_ERROR_56, but we no longer use the cable density to pre-check
    //         for non-convergence. Redo this error checking mechanism. 
    MAPSetUniversalErrorStat( code, str, err, msg );
  }
  
  SetCableDistance(); // Cable excursion in the element frame x and y directions
  InitializeHAndV( g , rho_sea ); // set H and V (as well as the anchor/fairlead force) 
};


// ============================================================================
// InitializeCableProperty
//
// @input : g            --
// @input : rho          -- 
// ============================================================================
void Element::
InitializeCableProperty( const double &g       , 
                         const double &rho_sea ) 
{
  // cross-section area and weight per length.
  // @note : weight per length is converted from N/m to kN/m to obtain the correct units
  //         to bring it consistent with the fairlead forces and EA (axial line stiffness)
  this->A = MAP_CONST::PI*pow( (this->line_property->Diam.value/2) , 2 );  
  this->omega = 1e-3*( g*(this->line_property->MassDenInAir.value) - g*A*rho_sea ); 
  
  // Make sure the cable is not neutrally buoyant.
  // If the cable sensity matches the sea density, then the program will raise an error message.
  if ( fabs( omega ) <= 1e-3 ){
    throw MAP_ERROR_56;
  }

  // if OMIT_CONTACT and cable density is less than seawater, throw a warning.
  if( omega < 0.0 && OMIT_CONTACT_flag == true ) {
    throw MAP_WARNING_1;
  }
};


// ============================================================================
// ResetNodes
//
// set the variables 'sum_FX', 'sum_FY' and 'sum_FZ' for the 'fairlead' and 
// 'anchor' nodes to zero. This method calls Node::SetSumForceToZero()
// ============================================================================
void Element::
ResetNodes( )
{
  this->fairlead->SetSumForceToZero( );    // for fairlead
  this->anchor->SetSumForceToZero  ( );    // for anchor
};


// ============================================================================
// UpdateElement
//
// The element is updated at each step in the solve process. We update:
//   -- psi
//   -- the cable space (l) and height (h)
//   -- also, as psi, l and h are updated, the fairlead/anchor forces in the 
//      element frame changes. This must be updated as well. 
// @todo : remove in future versions once we are satisfied this behaves correctly
// ============================================================================
void Element::
UpdateElement( MAP_ErrStat_class &err , 
               MAP_Message_class &msg ) 
{ 
  this->SetCableDistance( ); 
  
  // set the angle between global and local element axis
  this->SetPsi( );

  if ( this->fairlead->FX.is_fixed==true && this->fairlead->FY.is_fixed==false ){
    // if FX is fixed and FY is an iterated parameter:
    this->Hx = this->fairlead->FX.value;             
    this->Hy = tan( psi ) * this->fairlead->FX.value;
    this->H.value = sqrt( pow(this->Hx,2) + pow(this->Hy,2) );
  } else if ( this->fairlead->FX.is_fixed==false && this->fairlead->FY.is_fixed==true ) {
    // if FX is an iterated parameter and FY is fixed:
    this->Hy = this->fairlead->FY.value;
    this->Hx = this->fairlead->FY.value/tan( psi );
    
    // make sure the demoniator in the above statement is not small (small divisor error) 
    assert( this->psi >= 1e-3 ); // @todo
    
    this->H.value = sqrt( pow(this->Hx,2) + pow(this->Hy,2) );
  } else if ( this->fairlead->FX.is_fixed==true && this->fairlead->FY.is_fixed==true && this->fairlead->type != Connect ) {
    // if FX and FY are fixed and the fairlead is not a Connect node:
    this->Hx = this->fairlead->FX.value;
    this->Hy = this->fairlead->FY.value;
  } else { // the default case where both FX and FY are iterated
    this->Hx = this->H.value * cos( psi );
    this->Hy = this->H.value * sin( psi );   
  };// END if

  this->Vz = this->V.value; // since z and Z are in alignment, this is true at all times 

  // set anchor forces
  if ( this->is_resting_on_seabed == false ) { // element is not resting on the seabed
    this->Ha_x = this->Hx;
    this->Ha_y = this->Hy;
    this->Va = Vz - omega*Lu.value;
  } else { // element is resting on the seabed
    // @todo : This algorithm is really, really ugly. Fix this.
    this->Lb = this->Lu.value - this->V.value/this->omega;

    this->Ha_x = (this->H.value - Lb*this->GetCB()*this->omega ) > 0 ? (this->H.value - Lb*this->GetCB()*this->omega )*cos( psi ) : 0;
    this->Ha_y = (this->H.value - Lb*this->GetCB()*this->omega ) > 0 ? (this->H.value - Lb*this->GetCB()*this->omega )*sin( psi ) : 0;
    this->Va = 0.0;
  };// END if

  this->HX.value = this->Hx; // used solely to print to output file:
  this->HY.value = this->Hy; // used solely to print to output file:

  this->AddForceToFairleadNode(); // add the fairlead forces to the node
  this->AddForceToAnchorNode();   // add the fairlead forces to the node
};


// ============================================================================
// SetPsi
// 
// 'psi' is the angle between the element x-axis (local frame) and X-axis 
// (global frame). This essentially produces this rotation matrix:
//
//           | cos(psi)  -sin(psi)  0 |
// R(psi) =  | sin(psi)   cos(psi)  0 |
//           |    0         0       1 |
//
//      1) first find psi - the angle of rotation between the element frame and the
//      2) global reference frame
//      3) r_j = fairlead displacement
//      4) r_i = anchor displacement
//      5) acos( dot( (r_j-r_i) , (u_i) ) / ( norm(r_j-r_i) ) )
// ============================================================================
void Element::
SetPsi( ) 
{
  // check if the anchor and fairlead nodes occupy the same position in space
  double overlap = sqrt( pow( (this->fairlead->X.value - this->anchor->X.value) , 2) + 
                         pow( (this->fairlead->Y.value - this->anchor->Y.value) , 2) + 
                         pow( (this->fairlead->Z.value - this->anchor->Z.value) , 2) );
  
  // make sure the demoninator is not zero (or close to it) 
  if ( overlap <= 1e-2 ) { throw MAP_ERROR_55; };

  this->norm = sqrt( pow( (this->fairlead->X.value - this->anchor->X.value) , 2) + 
                     pow( (this->fairlead->Y.value - this->anchor->Y.value) , 2) );

//  // This is done for the exceptional case of a perfectly vertical cable
//  // @note : if you have problems dealing with a perfectly vertical cable, 
//  //         this is the culprit.
//  if( this->norm<=1e-2) { 
//    this->norm = 1e-2; 
//  } 

  // find the angle Psi
  if ( (this->fairlead->Y.value - this->anchor->Y.value) >= 0) {
    // this simply finds the angle psi between the local and
    // global reference frames simply by evaluation trig relationships
    this->psi = acos( (this->fairlead->X.value - this->anchor->X.value)/this->norm );
  } else {
    this->psi = -acos( (this->fairlead->X.value - this->anchor->X.value)/this->norm );
  }
};


// ============================================================================
// SetCableDistance
//
// l -- the cable span is equal to the norm of the difference between the 
//      fairlead and anchor displacement vectors (x and y components only) 
//
// h -- the z and Z axis are in alignment; so the cable height is equal to the 
//      difference between the anchor and fairlead (z component only)
// ============================================================================
inline void Element::
SetCableDistance( ) 
{
  // set l (horizontal displacement between fairlead and anchor) and h (vertical
  // displacement between fairlead and anchor)
  this->l  = sqrt( pow( (this->fairlead->X.value - this->anchor->X.value) , 2) + 
                   pow( (this->fairlead->Y.value - this->anchor->Y.value) , 2) );
  this->h = this->fairlead->Z.value - this->anchor->Z.value;

//  // This is done for the exceptional case of a perfectly vertical cable
//  // @note : if you have problems dealing with a perfectly vertical cable, 
//  //         this is the culprit.
//  if ( this->l < 1e-2 ) {
//    this->l = 1e-2;
//    //this->H.value = 1e-2;
//  }
};


// ============================================================================
// InitializeHAndV
//
// If H and V are not user supplied initial guesses or fixed values, then the 
// process to estimate H and V is invoked using the method given in:
//
// @article{peyrot1979,
//     title={{Analysis of Cable Structures}},
//     author={Peyrot, A. H., and Goulois, A. M.},
//     journal={Computers and Structures},
//     volume={10},
//     pages={805--813},
//     year={1979},
//     publisher={Pergamon Press, Ltd.}
// }
//
// @input : g            --
// @input : rho          --
// ============================================================================
void Element::
InitializeHAndV( const double &g       , 
                 const double &rho_sea ) 
{
  // initialize lambda
  double lambda = 9999.9;
  double TOL = 1e-2;

  norm = sqrt( pow( (fairlead->X.value - anchor->X.value) , 2) + 
               pow( (fairlead->Y.value - anchor->Y.value) , 2) );
  
  // find the angle psi
  if ( norm <= TOL ) {
    psi = 1.0;//checkpoint();
  } else {
    if ( (fairlead->Y.value - anchor->Y.value) >= 0) {
      // this simply finds the angle psi between the local and
      // global reference frames simply by evaluation trig relationships
      psi = acos( (fairlead->X.value - anchor->X.value)/norm );
    } else {
      psi = -acos( (fairlead->X.value - anchor->X.value)/norm );
    }
  }
  
  if ( l==0 ){
    lambda = 1000000;
  } else if ( sqrt( pow(l,2) + pow(h,2) ) >= this->Lu.value ) {
    lambda = 0.2;
  } else {
    lambda = sqrt(3*( (pow(this->Lu.value,2) - pow(h,2))/pow(l,2) - 1));
  }

  assert( lambda != 9999.9 );

  // set up temporary variables to store the estimate forces from Peyrot & Goulois
  const double tempH = fabs( omega*l / ( 2*lambda ) );
  const double tempV = ( omega/2 ) * ( (h/tanh( lambda ) ) + Lu.value );

  // fairlead node forces
  if ( this->fairlead->type != Connect ) {

    // =======  Set X direction fairlead force  ====== 
    if ( this->fairlead->FX.is_fixed==true ){
      Hx = this->fairlead->FX.value;
    } else { 
      Hx = tempH * cos( this->psi );
    }
    
    // =======  Set Y direction fairlead force  ====== 
    if ( this->fairlead->FY.is_fixed==true ) {
      Hy = this->fairlead->FY.value;
    } else {
      Hy = tempH * sin( this->psi );
    }
    
    // =======  Set Z direction fairlead force  ====== 
    if ( this->fairlead->FZ.is_fixed==true ) { // if FZ has a fixed value 
      Vz = this->fairlead->FZ.value - this->fairlead->M.value*g 
        + this->fairlead->B.value*( g*rho_sea );
    } else {
      Vz = tempV;
    }
  } else { // this is the default for a connect node    
    Hx = tempH * cos( this->psi ); // This is H in the local element frame
    Hy = tempH * sin( this->psi ); // This is H in the local element frame
    Vz = tempV;                    // This is V in the local element frame
  }

  // the value for H (in the element frame) is equal to the vector 
  // sum of Hx and Hy
  this->V.value = Vz;
  this->H.value = sqrt( pow(Hx,2) + pow(Hy,2) );

  if( H.value <= TOL ) H.value = TOL;

  // anchor node forces
  if (this->anchor->type != Connect) {
    // =======  Set X direction anchor force  ======
    if ( this->anchor->FX.is_fixed==true ) { // if the FX has a fixed value
      Ha_x = this->anchor->FX.value;
    } else { 
      Ha_x = H.value * cos( this->psi );
    }
        
    // =======  Set Y direction anchor force  ====== 
    if ( this->anchor->FY.is_fixed==true ) { // if the FY has a fixed value
      Ha_y = this->anchor->FY.value;
    } else {
      Ha_y = H.value * sin( this->psi );
    }

    // =======  Set Z direction anchor force  ====== 
    // if the FZ has a fixed value
    if ( this->anchor->FZ.is_fixed==true ) {
      Va = this->anchor->FZ.value;
      Va = this->anchor->FZ.value - this->anchor->M.value*g
        + this->anchor->B.value*( g*rho_sea );
    } else {
      Va = Vz - omega*Lu.value;
    }
  } else { // this is the default for a connect node        
    Ha_x = H.value * cos( this->psi ); // This is H in the local element frame
    Ha_y = H.value * sin( this->psi ); // This is H in the local element frame
    Va = Vz - omega*Lu.value;          // This is V in the local element frame
  }

  // ===============  Update Node Force  =============================
  // 
  // The last step is to update the FX, FY and FZ. Fi is the force 
  // applied to the node to hold it in static equilibrium. This is 
  // essentially Hx, Hy and V (minus the effects of the buoyancy 
  // module and added node weight)
  // 
  // ======  !!!!!!  ====== 
  // 
  // THE NEXT LINES CHANGE VALUES SEEN IN MAP_OtherStateType_class
  // 
  // ======  !!!!!!  ====== 
  // 
  // An anchor node cannot be attached be attached to a vessel; hence 
  // there is no Vessel conditional check like the fairlead node  	   
  // 
  // =================================================================
  this->AddForceToFairleadNode();
  this->AddForceToAnchorNode();
};


// ============================================================================
// AddForceToFairleadNode
// ============================================================================
inline void Element::
AddForceToFairleadNode( )
{  
  // now sum this force in with the other applied node forces
  this->fairlead->AddToSumFX( Hx );
  this->fairlead->AddToSumFY( Hy );
  this->fairlead->AddToSumFZ( Vz );
};


// ============================================================================
// AddForceToAnchorNode
//
// We add the negative of the element anchor force to the node (Node::sum_FX,
// Node::sum_FY, Node::sum_FY)
// ============================================================================
inline void Element::
AddForceToAnchorNode( )
{
  // the anchor foces are negative because they act opposite of the 
  // force applied at the fairlead
  this->anchor->AddToSumFX( -Ha_x );
  this->anchor->AddToSumFY( -Ha_y );
  this->anchor->AddToSumFZ( -Va );
};


// ============================================================================
// GetElementPlotArray
//
// Allows the cable to be plotted in Python using the matplotlib library
//
// @input  : X           -- X global position of element point [meters]
// @input  : Y           -- Y global position of element point [meters]
// @input  : Z           -- Z global position of element point [meters]
// @input  : Msg         -- Error message status
//
// @output : integer     -- error code. If '1', no error. If '-1', an error is 
//                          reported. 
// ============================================================================
bool Element::
GetElementPlotArray( std::string &X   , 
                     std::string &Y   ,
                     std::string &Z   ,
                     MAP_ErrStat_class &err ,
                     MAP_Message_class &msg ) 
{
  int    N     = 100;                              // number of points being plotted
  double n     = boost::lexical_cast<double>( N ); // number of points being plotted 
  double Fh    = this->H.value;                    // horizontal fairlead force
  double Fv    = this->V.value;                    // vertical fairlead force
  double w     = this->omega;			     // mass per length 
  double lu    = this->Lu.value;		     // unstretched length
  double ea    = this->GetEA();                    // element stiffness
  double angle = this->psi;                        // angles between element and global axis
  double cb    = this->GetCB();                    // cable/sea bed friction coefficient

  double x = 9999.9;                               // X global position of point on element
  double y = 9999.9;                               // Y global position of point on element
  double z = 9999.9;                               // Z global position of point on element
  double S = 9999.9;                               // position along cable element
    
  std::string temp_X = "";     
  std::string temp_Y = "";      
  std::string temp_Z = "";  
  std::string str = "";

  try {
    for ( int s=0 ; s<n+1 ; s++ ){
      S = (s/n)*lu;

      // If the cable is not resting on the seabed, we use the classic catenary equation
      // for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
      // the modified version as done in the FAST wind turbine program. 
      //
      // @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
      //        Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.

      // =======  Catenary equations for a hanging chain  ====== 
      if (this->is_resting_on_seabed == false ){ // Element is NOT resting on sea floor       
        x = ((Fh/w)*boost::math::asinh( Fv/Fh )   
             - (Fh/w)*boost::math::asinh( (Fv-S*w)/Fh )
             + (Fh*S)/(ea))*cos( angle ) - this->fairlead->X.value; // X position of element in global coordinates
      
        y = (  (Fh/w)*boost::math::asinh( Fv/Fh )     
               - (Fh/w)*boost::math::asinh( (Fv-S*w)/Fh )
               + (Fh*S)/(ea))*sin( angle ) - this->fairlead->Y.value; // Y position of element in global coordinates
      
        z =  (Fh/w)*( sqrt( 1+pow(Fv/Fh,2) ) - sqrt( 1+pow((Fv-w*S)/Fh,2) ) )
          + (1/(ea))*(Fv*S+w*S*S/2) - this->fairlead->Z.value; // Z position of element in global coordinates 
      } else { // Element is resting on sea floor
        Lb            = lu - (Fv/w);
        double lambda = 0.0; 

        // Find the right equation for X(s) and Y(s) to use for plotting the element         
        // resting on the seabed                                                             
        if ( 0<=S && S<=(Lb-Fh/(cb*w)) ) { // for 0 <= s <= Lb - H/(Cb*w)                            
          x = -(S) * cos( angle ) - this->anchor->X.value; // X position of element in global coordinates        
          y = -(S) * sin( angle ) - this->anchor->Y.value; // Y position of element in global coordinates
        } else if ( (Lb-Fh/(cb/w))<S && S<=Lb )  { // for Lb - H/(Cb*w) < s <= Lb              
          lambda = (Lb-Fh/(cb*w))>0 ? (Lb-Fh/(cb*w)) : 0;                                    
        
          x = -( S + ((cb*w)/(2*ea)) * (S*S - 2*(Lb-Fh/(cb*w))*S + (Lb- Fh/(cb*w))*lambda ) )
            * cos( angle ) - this->anchor->X.value; // X position of element in global coordinates 

          y = -( S + ((cb*w)/(2*ea)) * (S*S - 2*(Lb-Fh/(cb*w))*S + (Lb- Fh/(cb*w))*lambda ) )
            * sin( angle ) - this->anchor->Y.value; // Y position of element in global coordinates 
        } else { // for Lb < s <= L
          lambda = (Lb-Fh/(cb*w))>0 ? (Lb-Fh/(cb*w)) : 0;
        
          x = -( Lb + (Fh/w)*boost::math::asinh( (w*(S-Lb))/Fh ) ) * cos( angle )    
            + ( ((Fh*S)/(ea)) + ((cb*w)/(2*ea))*( -Lb*Lb + (Lb-Fh/(cb*w))*lambda ) ) 
            * cos( angle ) - this->anchor->X.value; // X position of element in global coordinates 

          y = -( Lb + (Fh/w)*boost::math::asinh( (w*(S-Lb))/Fh ) ) * sin( angle )
            + ( ((Fh*S)/(ea)) + ((cb*w)/(2*ea))*( -Lb*Lb + (Lb-Fh/(cb*w))*lambda ) )
            * sin( angle ) - this->anchor->Y.value;  // Y position of element in global coordinates  
        };

        // Find the right equation for Z(s) to use for plotting the element resting on the seabed 
        if ( 0<=S && S<=Lb ) {           
          z = 0 - this->anchor->Z.value; // Z position of element in global coordinates
        } else {        
          z = -( (Fh/w)*( sqrt(1 + pow((w*(S-Lb)/Fh),2) ) - 1 ) + ((w*pow((S-Lb),2))/(2*ea)) )
            - this->anchor->Z.value; // Z position of element in global coordinates
        };
      };

      // Now check to make sure the valid numbers are calculated for the x, y and z 
      // direction. The following lines of code convert real numbers into string arguments
      // for plotting purposes. 
      try { // X direction	    
        X += boost::lexical_cast<std::string>( -x );
        if ( s != n ) X += " , ";
      } catch ( boost::bad_lexical_cast const& ) {
        MAPSetUniversalErrorStat( MAP_ERROR_78 , str, err, msg );
        return false;
      }

      try { // Y direction
        Y += boost::lexical_cast<std::string>( -y );
        if ( s != n ) Y += " , ";
      } catch ( boost::bad_lexical_cast const& ) {
        MAPSetUniversalErrorStat( MAP_ERROR_79 , str, err, msg );
        return false;
      }

      try { // Z direction
        Z += boost::lexical_cast<std::string>( -z );
        if ( s != n ) Z += " , ";
      } catch ( boost::bad_lexical_cast const& ) {
        MAPSetUniversalErrorStat( MAP_ERROR_80 , str, err, msg );
        return false;
      }
    }
  } catch( std::exception const& excep ){ 
    std::string errStr = "";
    std::ostringstream S;
    S << this->H.index;
    errStr = S.str();
    MAPSetUniversalErrorStat( MAP_ERROR_87 , errStr, err, msg );    
    return false;
  }
  return true;
};


/**
 *
 *
 */
void Element::
GetElementPlotArray( std::vector<double> &arrX , 
                     std::vector<double> &arrY ,
                     std::vector<double> &arrZ ,
                     MAP_ErrStat_class   &err  ,
                     MAP_Message_class   &msg  ) 
{
  int    N     = 100;                              // number of points being plotted
  double n     = boost::lexical_cast<double>( N ); // number of points being plotted 
  double Fh    = this->H.value;                    // horizontal fairlead force
  double Fv    = this->V.value;                    // vertical fairlead force
  double w     = this->omega;			   // mass per length 
  double lu    = this->Lu.value;		   // unstretched length
  double ea    = this->GetEA();                    // element stiffness
  double angle = this->psi;                        // angles between element and global axis
  double cb    = this->GetCB();                    // cable/sea bed friction coefficient

  double x = 9999.9;                               // X global position of point on element
  double y = 9999.9;                               // Y global position of point on element
  double z = 9999.9;                               // Z global position of point on element
  double S = 9999.9;                               // position along cable element

  std::cout << this->fairlead->Z.value << std::endl;    
  for ( int s=0 ; s<n+1 ; s++ ){
    S = (s/n)*lu;

    /**
     * If the cable is not resting on the seabed, we use the classic catenary equation
     * for a hanging chain to plot the mooring line profile. Otherwise if it is, we 
     * the modified version as done in the FAST wind turbine program. 
     *
     * @ref : J. Jonkman, November 2007. "Dynamic Modeling and Loads Analysis of an 
     *        Offshore Floating Wind Turbine." NREL Technical Report NREL/TP-500-41958.
     */

    // Catenary equations for a hanging chain
    if (this->is_resting_on_seabed == false ){ 
      // Element is NOT resting on sea floor       
      // X/Y position of element in global coordinates
      x = ((Fh/w)*boost::math::asinh( Fv/Fh ) - (Fh/w)*boost::math::asinh( (Fv-S*w)/Fh ) + (Fh*S)/(ea))*cos( angle ) - this->fairlead->X.value; 
      y = (  (Fh/w)*boost::math::asinh( Fv/Fh ) - (Fh/w)*boost::math::asinh( (Fv-S*w)/Fh ) + (Fh*S)/(ea))*sin( angle ) - this->fairlead->Y.value; 
      
      // Z position of element in global coordinates 
      z =  (Fh/w)*( sqrt( 1+pow(Fv/Fh,2) ) - sqrt( 1+pow((Fv-w*S)/Fh,2) ) ) + (1/(ea))*(Fv*S+w*S*S/2) - this->fairlead->Z.value; 
    } else { 
      // Element is resting on sea floor
      Lb            = lu - (Fv/w);
      double lambda = 0.0; 

      // Find the right equation for X(s) and Y(s) to use for plotting the element resting on the seabed                                                             
      if ( 0<=S && S<=(Lb-Fh/(cb*w)) ) { 
        // for 0 <= s <= Lb - H/(Cb*w)                            
        x = -(S) * cos( angle ) - this->anchor->X.value; // X position of element in global coordinates        
        y = -(S) * sin( angle ) - this->anchor->Y.value; // Y position of element in global coordinates
      } else if ( (Lb-Fh/(cb/w))<S && S<=Lb )  { 
        // for Lb - H/(Cb*w) < s <= Lb              
        lambda = (Lb-Fh/(cb*w))>0 ? (Lb-Fh/(cb*w)) : 0;                                    
        
        // X/Y position of element in global coordinates 
        x = -( S + ((cb*w)/(2*ea)) * (S*S - 2*(Lb-Fh/(cb*w))*S + (Lb- Fh/(cb*w))*lambda ) ) * cos( angle ) - this->anchor->X.value; 
        y = -( S + ((cb*w)/(2*ea)) * (S*S - 2*(Lb-Fh/(cb*w))*S + (Lb- Fh/(cb*w))*lambda ) ) * sin( angle ) - this->anchor->Y.value; 
      } else { 
        // for Lb < s <= L
        lambda = (Lb-Fh/(cb*w))>0 ? (Lb-Fh/(cb*w)) : 0;

        // X/Y position of element in global coordinates         
        x = -( Lb + (Fh/w)*boost::math::asinh( (w*(S-Lb))/Fh ) ) * cos( angle )    
          + ( ((Fh*S)/(ea)) + ((cb*w)/(2*ea))*( -Lb*Lb + (Lb-Fh/(cb*w))*lambda ) ) * cos( angle ) - this->anchor->X.value; 
        y = -( Lb + (Fh/w)*boost::math::asinh( (w*(S-Lb))/Fh ) ) * sin( angle )
          + ( ((Fh*S)/(ea)) + ((cb*w)/(2*ea))*( -Lb*Lb + (Lb-Fh/(cb*w))*lambda ) ) * sin( angle ) - this->anchor->Y.value;  
      }

      // Find the right equation for Z(s) to use for plotting the element resting on the seabed 
      if ( 0<=S && S<=Lb ) {           
        // Z position of element in global coordinates
        z = 0 - this->anchor->Z.value; 
      } else {        
        // Z position of element in global coordinates
        z = -( (Fh/w)*( sqrt(1 + pow((w*(S-Lb)/Fh),2) ) - 1 ) + ((w*pow((S-Lb),2))/(2*ea)) ) - this->anchor->Z.value; 
      }
    }
    arrX.push_back(-x);
    arrY.push_back(-y);
    arrZ.push_back(-z);
  }
};



// ============================================================================
// GetLineTensionString
//
// Prints the tension in the fairlead to the MAP map.out file
//
// @output : string           -- tension in for 10 points along the cable 
//                               [kilo-Newtons]
//
// @todo :                    -- right now, the code is hard coded to evaluate 
//                               10 points along the line. In the future, we 
//                               may want to change this so an arbitrary
//                               number of points are evaluated
// ============================================================================
std::string Element:: 
GetLineTensionString( ) 
{
  std::string        out         = ""; // initlize the string we are writting for the outputs
  std::string        temp        = "";
  std::ostringstream S;
  double             length      = 9999.9;
  double             temp_double = 9999.9;
    
  S << std::fixed << std::setprecision(3);
    
  double Fh    = this->H.value;                    // horizontal fairlead force
  double Fv    = this->V.value;                    // vertical fairlead force
  double w     = this->omega;			     // mass per length 
  double cb    = this->GetCB();                    // cable/sea bed friction coefficient
  double lu    = this->Lu.value;		     // unstretched length
  double Lb    = lu - (Fv/w);			     
    
  for (int i=0 ; i<10 ; i++){ // i<10 since we are only evaluating 10 points along cable

    // @todo : right now, the code is  hard coded to evaluate 10 points along the line. In the
    //         future, we may want to change this so an arbitrary number of points are evaluated
    length = (static_cast<double>(i)/9)*(this->Lu.value); 

    // =======  Cases for cable hanging in air or touching seabed  ======     <----------------------------------+
    if ( this->is_resting_on_seabed==false  ) { // Element is NOT resting on sea floor
      S << sqrt( pow((this->H.value),2) + pow((this->Va+this->omega*length),2) );
      temp += S.str();
      S.str("");S.clear();

      // make sure the string array is long enough to put at least one white space
      // at the tail; otherwise, the do-while loop will result in a program crash
      assert( temp.size() < 14 );

      // fill string with white spaces to bring it to the correct length
      do { temp += " "; } while( temp.size()<15 ); // @todo : 15 is arbitrarily chosen.

      out += temp;
      temp = "";  
    } else { // cable is touching the seabed 
      if ( 0<=length && length<=Lb ){ // for 0 <= s <= Lb
        temp_double = (Fh+cb*w*(length-Lb))>0 ? (Fh+cb*w*(length-Lb)) : 0;
        S << temp_double;
      } else { // for Lb < s <= L
        S << sqrt( Fh*Fh + pow( (w*(length-Lb)) ,2) );
      };

      temp += S.str();    
      S.str("");S.clear();

      // make sure the string array is long enough to put at least one white space
      // at the tail; otherwise, the do-while loop will result in a program crash
      assert( temp.size() < 14 );                         

      // fill string with white spaces to bring it to the correct length
      do { temp += " "; } while( temp.size()<15 ); // @todo : 15 is arbitrarily chosen. 

      out += temp;           
      temp = "";
    }
  }

  return out;
};


// ============================================================================
// GetLineTensionString
//
// Prints the tension in the fairlead to the MAP map.out file
//
// @output : string           -- tension in for 10 points along the cable 
//                               [Newtons]
//
// @todo :                    -- right now, the code is hard coded to evaluate 
//                               10 points along the line. In the future, we 
//                               may want to change this so an arbitrary
//                               number of points are evaluated
// ============================================================================
void Element:: 
GetLineTensionValues( float * passArr ) 
{
  float length      = 999.9;
  float temp_double = 999.9;    
  float Fh          = static_cast<float>( this->H.value  );   // horizontal fairlead force
  float Fv          = static_cast<float>( this->V.value  );   // vertical fairlead force
  float w           = static_cast<float>( this->omega    );   // mass per length 
  float cb          = static_cast<float>( this->GetCB()  );   // cable/sea bed friction coefficient
  float lu          = static_cast<float>( this->Lu.value );   // unstretched length
  float Lb          = static_cast<float>( lu - (Fv/w)    );
    
  for (int i=0 ; i<10 ; i++){ // i<10 since we are only evaluating 10 points along cable
    
    // @todo : right now, the code is  hard coded to evaluate 10 points along the line. In the
    //         future, we may want to change this so an arbitrary number of points are evaluated
    length = (static_cast<float>(i)/9)*(this->Lu.value); 

    // =======  Cases for cable hanging in air or touching seabed  ======     <----------------------------------+
    if ( this->is_resting_on_seabed==false  ) { // Element is NOT resting on sea floor
      passArr[i] = sqrt( pow((this->H.value),2) + pow((this->Va+this->omega*length),2) );
    } else { // cable is touching the seabed 
      if ( 0<=length && length<=Lb ){ // for 0 <= s <= Lb
        passArr[i] = (Fh+cb*w*(length-Lb))>0 ? (Fh+cb*w*(length-Lb)) : 0;        
      } else { // for Lb < s <= L
        passArr[i] =  sqrt( Fh*Fh + pow( (w*(length-Lb)) ,2) );
      }
    }
  }
};


// ============================================================================
// GetLineTensionStringHeader 
//
// Prints the heades "T[1] T[2] ... T[n]" header for the element in the MAP
// map.out file
//
// @input  : int i             -- the element number
//
// @output : string            -- "T[1] T[2] ... T[n]"
// ============================================================================
std::string Element::
GetLineTensionStringHeader( const int i ) 
{
  std::string        out  = ""; // initialize the string we are writing for the outputs
  std::string        temp = "";
  std::ostringstream S;
  int NINE = 9;

  S << std::fixed << std::setprecision(1);

  for (int j=0 ; j<10 ; j++){
    // create the string we are writting as an output
    S << i+1;
    temp += "T[" + S.str() + "](";         
    S.str("");S.clear();

    // @todo : right now, the code is hard coded to evaluate 10 points along the line.
    // In the future, we may want to change this so an arbitrary number of points are 
    // evaluated
    S <<  ( static_cast<double>(j)/NINE )*(this->Lu.value);
 
    temp += S.str();
    temp += ")";
    S.str("");S.clear();

    // make sure the string array is long enough to put at least one white space 
    // at the tail; otherwise, the do-while loop will result in a program crash
    assert( temp.size() < 14 );

    // fill string with white spaces to bring it to the correct length
    do { temp += " "; } while( temp.size()<15 ); // @note : 15 is arbitrarily chosen. Fix this ? 

    out += temp;
    temp = "";
  };// END for 
    
  return out;
};


// ============================================================================
// GetLineTensionStringHeader 
//
// Prints the heades "T[1] T[2] ... T[n]" header for the element in the MAP
// map.out file
//
// @input  : int i             -- the element number
//
// @output : string            -- "T[1] T[2] ... T[n]"
// ============================================================================
std::string Element::
GetLineTensionStringHeaderForFast( const int i, const int j ) 
{
  std::string        out  = ""; // initialize the string we are writing for the outputs
  std::string        temp = "";
  std::ostringstream S;
  int NINE = 9;

  S << std::fixed << std::setprecision(1);

  // create the string we are writting as an output
  S << i+1;
  temp += "T[" + S.str() + "](";         
  S.str("");S.clear();
  
  // @todo : right now, the code is hard coded to evaluate 10 points along the line.
  // In the future, we may want to change this so an arbitrary number of points are 
  // evaluated
  S <<  ( static_cast<double>(j)/NINE )*(this->Lu.value);
  
  temp += S.str();
  temp += ")";
  S.str("");S.clear();
  
  // make sure the string array is long enough to put at least one white space 
  // at the tail; otherwise, the do-while loop will result in a program crash
  assert( temp.size() < 14 );
  
  // fill string with white spaces to bring it to the correct length
  do { temp += " "; } while( temp.size()<15 ); // @note : 15 is arbitrarily chosen. Fix this ? 
  
  out += temp;
      
  return out;
};


// ============================================================================
// GetLayLengthStringHeader 
// ============================================================================
std::string Element::
GetLayLengthStringHeader( const int i ) 
{
  std::string        out  = ""; // initialize the string we are writing for the outputs
  std::string        temp = "";
  std::ostringstream S;
  int NINE = 9;

  S << std::fixed << std::setprecision(1);

  // create the string we are writting as an output
  S << i+1;
  temp += "Lb[" + S.str() + "]";         
  S.str("");S.clear();

  // fill string with white spaces to bring it to the correct length
  do { temp += " "; } while( temp.size()<15 ); // @note : 15 is arbitrarily chosen. Fix this ? 
  
  out += temp;
  temp = "";
    
  return out;
};


/**
 *
 */
std::string Element::
GetFairleadTensionStringHeader( const int i ) 
{
  std::string        out  = ""; // initialize the string we are writing for the outputs
  std::string        temp = "";
  std::ostringstream S;
  int NINE = 9;

  S << std::fixed << std::setprecision(1);

  // create the string we are writting as an output
  S << i+1;
  temp += "TFair[" + S.str() + "]";         
  S.str("");S.clear();

  // fill string with white spaces to bring it to the correct length
  do { temp += " "; } while( temp.size()<15 ); // @note : 15 is arbitrarily chosen. Fix this ? 
  
  out += temp;
  temp = "";
    
  return out;
};


/**
 *
 */
std::string Element::
GetAnchorTensionStringHeader( const int i ) 
{
  std::string        out  = ""; // initialize the string we are writing for the outputs
  std::string        temp = "";
  std::ostringstream S;
  int NINE = 9;

  S << std::fixed << std::setprecision(1);

  // create the string we are writting as an output
  S << i+1;
  temp += "TAnch[" + S.str() + "]";         
  S.str("");S.clear();

  // fill string with white spaces to bring it to the correct length
  do { temp += " "; } while( temp.size()<15 ); // @note : 15 is arbitrarily chosen. Fix this ? 
  
  out += temp;
  temp = "";
    
  return out;
};


/**
 *
 */
double Element::
GetFairleadTensionMagnitude( ) const
{
  return sqrt( Hx*Hx + Hy*Hy + Vz*Vz );
}


/**
 *S << std::fixed << std::setprecision(3);
S << sqrt( pow((this->H.value),2) + pow((this->Va+this->omega*length),2) );
S.str("");S.clear();
 */
double Element::
GetAnchorTensionMagnitude( ) const
{
  return sqrt( Ha_x*Ha_x + Ha_y*Ha_y + Va*Va );
};


/**
 *
 */
std::string Element::
GetFairleadTensionMagnitudeString( ) const
{
  std::ostringstream S;
  std::string out = "";
  S << std::fixed << std::setprecision(3);
  S << sqrt( Hx*Hx + Hy*Hy + Vz*Vz );
  out += S.str();
  do { out += " "; } while( out.size()<15 );
  return out;
}


/**
 *
 */
std::string Element::
GetAnchorTensionMagnitudeString( ) const
{
  std::ostringstream S;
  std::string out = "";
  S << std::fixed << std::setprecision(3);
  S << sqrt( Ha_x*Ha_x + Ha_y*Ha_y + Va*Va );
  out += S.str();
  do { out += " "; } while( out.size()<15 );
  return out;
};


// ============================================================================
// GetLayLengthString
// ============================================================================
std::string Element:: 
GetLayLengthString( ) 
{
  std::string        out         = ""; // initlize the string we are writting for the outputs
  std::ostringstream S;
  this->Lb = (this->Lu.value - this->Vz/this->omega);

  if ( this->Lb < 0.0 ) {
    this->Lb = 0.0;
  }

  S << this->Lb;
  out += S.str();
  S.str("");S.clear();

  // make sure the string array is long enough to put at least one white space
  // at the tail; otherwise, the do-while loop will result in a program crash
  assert( out.size() < 14 );
  
  // fill string with white spaces to bring it to the correct length
  do { out += " "; } while( out.size()<15 ); // @todo : 15 is arbitrarily chosen.

  return out;
};


// ============================================================================
// GetLayLengthString
// ============================================================================
float Element:: 
GetLayLengthValue( ) const
{
  if ( this->Lb < 0.0 ) {
    return 0.0 ;
  } else {
    return static_cast<float>(this->Lb);
  }
};


// ============================================================================
// CompareNodeAddressWithFairlead
// ============================================================================
bool Element::
CompareNodeAddressWithFairlead( const Node &ref ) const 
{ 
  if ( &ref==this->fairlead ) {
    return true;
  } else {
    return false;
  }
};


// ============================================================================
// CompareNodeAddressWithAnchor
// ============================================================================
bool Element::
CompareNodeAddressWithAnchor( const Node &ref ) const 
{    
  if ( &ref==this->anchor ) {
    return true;
  } else {
    return false;
  };
};


// ============================================================================
// 
// ============================================================================
double Element::
GetPsi( ) const 
{ 
  return this->psi; 
}


// ============================================================================
// 
// ============================================================================
double Element::
GetOmega( ) const 
{ 
  return this->omega;                   
};


// ============================================================================
// 
// ============================================================================
double Element::
GetHeight( ) const 
{ 
  return this->h;                       
};


// ============================================================================
// 
// ============================================================================
double Element::
GetLength( ) const 
{ 
  return this->l;                       
};


// ============================================================================
// 
// ============================================================================
double Element::
GetArea( ) const 
{ 
  return this->A;                       
};


// ============================================================================
// 
// ============================================================================
double Element::
GetEA( ) const
{
  return this->line_property->EA.value; 
};


// ============================================================================
// 
// ============================================================================
double Element::
GetCB( ) const
{ 
  return this->line_property->CB.value; 
};


// ============================================================================
// 
// ============================================================================
bool Element::
GetRestingSeabedFlag( ) const 
{
  return is_resting_on_seabed;          
};


// ============================================================================
// 
// ============================================================================
void Element::
SetRestingSeabedFlag( const bool flag )
{
  this->is_resting_on_seabed = flag;    
};


// ============================================================================
// 
// ============================================================================
double Element::
GetXf( ) const 
{
  return this->GetFairleadPostion( &Node::X ); 
};


// ============================================================================
// 
// ============================================================================
double Element::
GetXa( ) const 
{ 
  return this->GetAnchorPosition  ( &Node::X ); 
};


// ============================================================================
// 
// ============================================================================
double Element::
GetYf( ) const 
{
  return this->GetFairleadPostion( &Node::Y ); 
};


// ============================================================================
// 
// ============================================================================
double Element::
GetYa( ) const 
{ 
  return this->GetAnchorPosition  ( &Node::Y ); 
};


// ============================================================================
// 
// ============================================================================
void Element::
SetHFlagTo( const bool flag )
{
  this->H.is_fixed = flag; 
}


// ============================================================================
// 
// ============================================================================
void Element::
SetVFlagTo( const bool flag )
{
  this->V.is_fixed = flag; 
}


// ============================================================================
// 
// ============================================================================
void Element::
SetPlotFlag( const bool flag ) 
{ 
  this->PLOT_flag = flag; 
}


// ============================================================================
// 
// ============================================================================
void Element::
SetXPosFlag( const bool flag ) 
{
  this->X_POS_flag = flag; 
}


// ============================================================================
// 
// ============================================================================
void Element::
SetYPosFlag( const bool flag ) 
{ 
  this->Y_POS_flag = flag; 
}


// ============================================================================
// 
// ============================================================================
void Element::
SetZPosFlag( const bool flag ) 
{
  this->Z_POS_flag = flag;
}


// ============================================================================
// 
// ============================================================================
void Element::
SetXForceFlag( const bool flag )
{
  this->X_FORCE_flag = flag; 
}


// ============================================================================
// 
// ============================================================================
void Element::
SetyForceFlag( const bool flag )
{ 
  this->Y_FORCE_flag  = flag; 
}


// ============================================================================
// 
// ============================================================================
void Element::
SetZForceFlag( const bool flag ) 
{
  this->Z_FORCE_flag = flag;
}


// ============================================================================
// 
// ============================================================================
void Element::
SetLineTensionFlag( const bool flag )
{
  this->LINE_TENSION_flag = flag;
}


/**
 *
 */
void Element::
SetFairTensionFlag( const bool flag )
{
  this->FAIR_TENSION_flag = flag;
}


/**
 *
 */
void Element::
SetAnchTensionFlag( const bool flag )
{
  this->ANCH_TENSION_flag = flag;
}


// ============================================================================
// 
// ============================================================================
void Element::
SetOmitContactFlag( const bool flag )
{
  this->OMIT_CONTACT_flag = flag; 
}


// ============================================================================
// 
// ============================================================================
void Element::
SetLayLengthFlag( const bool flag )
{
  this->LAY_LENGTH_flag = flag; 
}


// ============================================================================
// 
// ============================================================================
bool Element::
GetPlotFlag( ) const 
{
  return this->PLOT_flag;     
}


// ============================================================================
// 
// ============================================================================
bool Element::
GetXPosFlag( ) const
{ 
  return this->X_POS_flag;        
}


// ============================================================================
// 
// ============================================================================
bool Element::
GetYPosFlag( ) const 
{
  return this->Y_POS_flag;
}


// ============================================================================
// 
// ============================================================================
bool Element::
GetZPosFlag( ) const 
{
  return this->Z_POS_flag;
}


// ============================================================================
// 
// ============================================================================
bool Element::
GetXForceFlag( ) const 
{
  return this->X_FORCE_flag;     
}


// ============================================================================
// 
// ============================================================================
bool Element::
GetYForceFlag( ) const 
{
  return this->Y_FORCE_flag; 
}


// ============================================================================
// 
// ============================================================================
bool Element::
GetZForceFlag( ) const 
{
  return this->Z_FORCE_flag; 
}


// ============================================================================
// 
// ============================================================================
bool Element::
GetLineTensionFlag( ) const 
{
  return this->LINE_TENSION_flag; 
}


// ============================================================================
// 
// ============================================================================
bool Element::
GetOmitContactFlag( ) const 
{
  return this->OMIT_CONTACT_flag; 
}


// ============================================================================
// 
// ============================================================================
bool Element::
GetLayLengthFlag( ) const 
{
  return this->LAY_LENGTH_flag;   
}


/**
 * Checks to make sure the element does not double back on itself. This is
 * a limitation of the closed for analytical solution. This is a fatal error
 * that cannot be resolved without changing the mooring properties in the MAP
 * input file. 
 * 
 * @change  This used to be called in Element::UpdateStates( ), but it was causing false errors
 *          as the varaibles were being iterated. 
 * @see     MAP_OtherStateType_cass::Solve( )  This member is called once every solve iteration.
 * @throws  MAP_ERROR_85  Double backing of the element is occuring (unstretched line length is too long)
 */
void Element::
CheckMaximumLineLength( ) 
{
  double Lmax = 0.0;

  if( this->omega > 0.0 ) { // true when the line will sink in fluid
    // Compute the maximum stretched length of the line with seabed interaction beyond which the line
    // would have to double-back on itself; here the line forms an "L" between the anchor and fairlead 
    // (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead)
    Lmax = this->l - (this->GetEA()/this->omega) + sqrt( (this->GetEA()/this->omega)*(this->GetEA()/this->omega) 
                                                         + 2.0*(this->h)*(this->GetEA()/this->omega) );

    if( this->Lu.value>=Lmax && this->GetOmitContactFlag()==false ) {  // true if the line is as long or longer than its
                                                                       // maximum possible value with seabed interaction        
//      throw MAP_ERROR_85;
    };
  }
}


// ============================================================================
// CheckIfCableIsVertical
//
// Checks to determine if the cable element is vertical. If it is, the 
// fairlead/anchor (non-connect node, if it exists) is displaced by epislon
// so that the problem is no longer ill-conitioned (we hope).
// ============================================================================
void Element::
CheckIfCableIsVertical()
{
  double TOL = 5e-2;
  double epsilon = 5e-2;
  double height = fabs( fairlead->Z.value - anchor->Z.value );
  double offset = sqrt( pow( (fairlead->X.value - anchor->X.value) , 2) + 
                        pow( (fairlead->Y.value - anchor->Y.value) , 2) );

  if ( offset <= TOL*height ) {
    is_vertical = true;//SetVerticalFlag( true );
    if ( this->fairlead->GetNodeType() != Connect ) {
      ghost_value_x = fairlead->X.value;
      ghost_value_y = fairlead->Y.value;
      
      fairlead->X.value += epsilon;
      fairlead->Y.value += epsilon;
      
      H.value = TOL*V.value;
    } else if ( anchor->GetNodeType() != Connect ) {
    } else {
      checkpoint();
      // need to throw error here
    }
  } else {
    is_vertical = false;//SetVerticalFlag( false );
  }
}


// ============================================================================
// ResetGhostProperties
//
// Sets the is_verical flag to true or false. This is the first step before
// the cable system is solved. Use to prevent a perfectly vertcal cable to 
//  be solved, which will lead to an ill-conditioned Jacobian.
//
// @input  :  flag, true or false, determining if the element is perfectly 
//            vertical
// @see    :  MAP_OtherStateTypr_class::Solve(...)
// ============================================================================
void Element::
ResetGhostProperties( )
{
  if (is_vertical == true ){
    fairlead->X.value = ghost_value_x;
    fairlead->Y.value = ghost_value_y;
  }
}
