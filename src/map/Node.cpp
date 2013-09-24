/**
 * ====================================================================================================
 *                              Node.cpp
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


#include "Node.h"


/**
 * Initializes the VarType according to what is given as an argument into this member
 *
 * @param  varStr  The variable string to identify it by a name. name is used when writting to outfile  
 * @param  type    Set the node type as 'connect', 'fix', or 'vessel'
 * @param  i       Set the node index number. This is sequential. 
 * @param  ptr     Node
 * @param  err     error code
 * @param  msg     error message
 */
void Node::
SetVarType( const std::string &varStr , 
            const std::string &type   , 
            int               i       , 
            VarType Node::*   ptr     , 
            MAP_ErrStat_class       &err    , 
            MAP_Message_class       &msg    ) 
{
  ((*this).*ptr).name = type; // set the name of the VarType 
  ((*this).*ptr).index = i+1; // set the index number of the VarType
  VarType::SetGenericVarType( ((*this).*ptr) , 
                              varStr         , 
                              err            , 
                              msg            ); 
}


/**
 * Initialize sum_FX, sum_FY, sum_FZ to 0.0
 */
void Node::
SetSumForceToZero( )
{
  this->sum_FX = 0.0;
  this->sum_FY = 0.0;
  this->sum_FZ = 0.0;
}


/**
 * This sets the initial conditions for the Fix and Vessel nodes. This adds 
 * sum_FX, sum_FY and sum_FZ into the node.  Note that a 'Connect' node is 
 * not a part of this process because it is assumed the connect node has known 
 * FX, FY and FZ applied forces
 */
void Node::
InitializeNode( )
{
  //If the fairlead if fixed or attached to a vessel, do this:
  if( this->GetNodeType()==Fix || this->GetNodeType()==Vessel ){
    this->FX.value = sum_FX;
    this->FY.value = sum_FY;
    this->FZ.value = sum_FZ + this->M.value*grav - this->B.value*grav*sea_density;
  }
}


/**
 * Sets the equilibrium displacement of each node after the first Solve call
 *
 * @see  MAP_OtherStateType_class::initializeEquilibriumNodePosition()
 */
void Node::
SetEquilibriumDisplacement( ) 
{
  if ( this->type == Vessel ){
    this->SetPriorXValue( );  
    this->SetPriorYValue( );  
    this->SetPriorZValue( );  
  }
}


/** 
 * Return the global force in the X direction [kN]
 *
 * @property  VarType  FX is a value that can be iterated or fixed (decided at run-time by the MAP input file)
 * @return    double   the force in the X direction, 
 */
double Node::
GetFX( ) const 
{
  return this->FX.value; 
};


/** 
 * Return the global force in the Y direction [kN]
 *
 * @property  VarType  FY is a value that can be iterated or fixed (decided at run-time by the MAP input file)
 * @return    double   the force in the Y direction, 
 */
double Node::
GetFY( ) const 
{
  return this->FY.value; 
};


/** 
 * Return the global force in the Z direction [kN]
 * 
 * @property  VarType  FZ is a value that can be iterated or fixed (decided at run-time by the MAP input file)
 * @return    double   the force in the Z direction, 
 */
double Node::
GetFZ( ) const 
{ 
  return this->FZ.value; 
};


/** 
 * Return the node mass [kg]
 *
 * @property  VarType  M is a value that can be iterated or fixed (decided at run-time by the MAP input file)
 * @return    double   the mass in this node
 */
double Node::
GetM( ) const 
{ 
  return this->M.value; 
}; 


/** 
 * Return the node volumetric displacement [m^3]
 *
 * @property  VarType  B is a value that can be iterated or fixed (decided at run-time by the MAP input file)
 * @return    double   the dispalcement of the node. Must be multiplied by rho*g to get the buoyant force
 */
double Node::
GetB( ) const 
{ 
  return this->B.value; 
};


/** 
 * @return  double  the gravitational constant [m/s^2]
 */
double Node::
GetGrav( ) const 
{ 
  return this->grav;
};


/** 
 * @return  double  the sea density [kg/m^3]
 */
double Node::
GetSeaDensity( ) const 
{
  return this->sea_density;
};


/** 
 * Get the value of a particular VarType. The VarType can be any one defined in the Node.h class.
 *
 * @see     Node::X 
 * @see     Node::Y 
 * @see     Node::Z 
 * @see     Node::M 
 * @see     Node::B
 * @see     Node::FZ 
 * @see     Node::FY 
 * @see     Node::FZ 
 * @param   VarType   ptr is a point to a Node Vartype. Can be either M, B, FX, FY, FZ, X, Y or Z.
 * @return  double    the value of the VarType. Units can be [m], [kN], [kg], or [m^3].
 */
double Node::
GetVarTypeValue( VarType Node::* ptr ) const   
{
  return (this->*ptr).value; 
}


/** 
 * @param   VarType  ptr is a point to a Node Vartype. Can be either M, B, FX, FY, FZ, X, Y or Z.
 * @return  double   the value of the VarType. Units can be [m], [kN], [kg], or [m^3].
 */
void Node::
SetVarTypeValue( VarType Node::* ptr   , 
                 double          value ) 
{ 
  (this->*ptr).value = value;    
};


/** 
 * The type of node (i.e., Connect, Vessel, Fix) is set in the Node constructor 
 *
 * @see     Node::Node()
 * @see     MAP_OtherStateType_class::$addNode
 * @return  NodeType  the node type. It can be a 'Fix', 'Connect', or 'Vessel'.
 */
NodeType Node::
GetNodeType( ) const 
{ 
  return this->type; 
};


/** 
 * Stores the current X.value before X.value is modified.
 *
 * @see  Node::SetEquilibriumPosition() 
 * @see  MAP_OtherStateType_call::initializeEquilibriumNodePosition()
 */
void Node::
SetPriorXValue( )       
{ 
  this->prior_value_x = this->X.value; 
};


/** 
 * Stores the current Y.value before Y.value is modified.
 *
 * @see  Node::SetEquilibriumPosition()
 * @see  MAP_OtherStateType_call::initializeEquilibriumNodePosition()
 */
void Node::
SetPriorYValue( ) 
{ 
  this->prior_value_y = this->Y.value; 
};


/** 
 * Stores the current Z.value before Z.value is modified.
 * 
 * @see  Node::SetEquilibriumPosition()
 * @see  MAP_OtherStateType_call::initializeEquilibriumNodePosition()
 */
void Node::
SetPriorZValue( )       
{ 
  this->prior_value_z = this->Z.value; 
};


/** 
 * Returns the previous X node displacement before being touched by the solver
 *
 * @see     MAP_OtherStateType_class::writeLinearizedStiffnessMatrix()
 * @return  double the X node position before it was modified [m]
 */
double Node::
GetPriorXValue( ) const 
{ 
  return this->prior_value_x;          
};


/** 
 * Returns the previous X node displacement before being touched by the solver
 *
 * @see     MAP_OtherStateType_class::writeLinearizedStiffnessMatrix()
 * @return  double the X node position before it was modified [m]
 */
double Node::
GetPriorYValue( ) const 
{ 
  return this->prior_value_y;          
};


/** 
 * Returns the previous X node displacement before being touched by the solver
 *
 * @see     MAP_OtherStateType_class::writeLinearizedStiffnessMatrix()
 * @return  double the X node position before it was modified [m]
 */
double Node::
GetPriorZValue( ) const 
{ 
  return this->prior_value_z;          
};
   

/**
 * Set the 'solve_X_Newton_equation' boolean. This idicates if we are solving the
 * Newton equation so that static equilibrium for this node. For Connect nodes,
 * the flag will be "true"
 *
 * @param  bool  flag
 * @todo   should solve_Z_Newton_equation be a const?
 */
void Node::
SetXNewtonEquationFlag( bool flag )
{ 
  this->solve_X_Newton_equation = flag; 
};


/**
 * Set the 'solve_Y_Newton_equation' boolean. This idicates if we are solving the
 * Newton equation so that static equilibrium for this node. For Connect nodes,
 * the flag will be "true"
 *
 * @param  bool  flag
 * @todo   should solve_Z_Newton_equation be a const?
 */
void Node::
SetYNewtonEquationFlag( bool flag )
{ 
  this->solve_Y_Newton_equation = flag; 
};


/**
 * Set the 'solve_Z_Newton_equation' boolean. This idicates if we are solving the
 * Newton equation so that static equilibrium for this node. For Connect nodes,
 * the flag will be "true"
 *
 * @param  bool  flag
 * @todo   should solve_Z_Newton_equation be a const?
 */
void Node::
SetZNewtonEquationFlag( bool flag ) 
{
  this->solve_Z_Newton_equation = flag; 
};


/**
 * Get the 'solve_X_Newton_equation' boolean. 
 *
 * @return  bool  flag indicating if we are solving the force-balance equation in X
 */
bool Node::
GetXNewtonEquationFlag( ) const 
{
  return this->solve_X_Newton_equation; 
};


/**
 * Get the 'solve_Y_Newton_equation' boolean. 
 *
 * @return  bool  flag indicating if we are solving the force-balance equation in Y
 */
bool Node::
GetYNewtonEquationFlag( ) const
{
  return this->solve_Y_Newton_equation; 
};


/**
 * Get the 'solve_Z_Newton_equation' boolean. 
 *
 * @return  bool  flag indicating if we are solving the force-balance equation in Z
 */
bool Node::
GetZNewtonEquationFlag( ) const 
{
  return this->solve_Z_Newton_equation; 
};


/**
 * Sum forces in the fairlead/anchor in each element to the node
 *
 * @see    Element::AddForceToFairleadNode( ) 
 * @see    Element::AddForceToAnchorNode( )
 * @param  double  $sum_FX the delta X force we incrementing the sum forces by
 */
void Node::
AddToSumFX( double x_force ) 
{ 
  this->sum_FX += x_force;
};


/**
 * Sum forces in the fairlead/anchor in each element to the node
 *
 * @see    Element::AddForceToFairleadNode( )
 * @see    Element::AddForceToAnchorNode( )
 * @param  double  $sum_FY the delta Y force we incrementing the sum forces by
 */
void Node::
AddToSumFY( double y_force )
{
  this->sum_FY += y_force;
};


/**
 * Sum forces in the fairlead/anchor in each element to the node
 *
 * @see    Element::AddForceToFairleadNode( ) 
 * @see    Element::AddForceToAnchorNode( )
 * @param  double  $sum_FZ the delta Z force we incrementing the sum forces by
 */
void Node::
AddToSumFZ( double z_force ) 
{ 
  this->sum_FZ += z_force;
};


/**
 * Pass the sum forces previous calculated in the node
 *
 * @see     SUM_F_IN_X() overloaded function
 * @see     Node::AddToSumFX()
 * @return  double  the node sum forces in X [kN]
 */
double Node::
GetSumFX( ) const 
{ 
  return this->sum_FX; 
};


/**
 * Pass the sum forces previous calculated in the node
 *
 * @see     SUM_F_IN_Y() overloaded function
 * @see     Node::AddToSumFY()
 * @return  double  the node sum forces in Y [kN]
 */
double Node::
GetSumFY( ) const 
{
  return this->sum_FY;
};


/**
 * Pass the sum forces previous calculated in the node
 *
 * @see     SUM_F_IN_Z() overloaded function
 * @see     Node::AddToSumFZ()
 * @return  double  the node sum forces in Z [kN]
 */
double Node::
GetSumFZ( ) const 
{
  return this->sum_FZ;  
};
