/**
 * ====================================================================================================
 *                              Node.h
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


#ifndef _NODE_H
#define _NODE_H


#include "VarType.h"   // Preprocessor Defitions in VarType.h
                       //   #include <boost/lexical_cast.hpp>
                       //   #include <boost/algorithm/string.hpp>
                       //   #include <string>
                       //   #include <iomanip>
                       //   #include "MAP_Message_class.h" 
                       //   #include "MAP_ErrStat_class.h" 

// Node predefition in order to create the Solving Function 
// structs. The Node object is declared further below.
struct Node;  


// ==========   Solving functions   ===============     <--------------------+
/** 
 * Newton equation (force balace) struct for X, Y and Z,  For each node, there 
 * are a maximum of three equations we are solving for:                             
 *   - X direction sum forces                   
 *   - Y direction sum forces                   
 *   - Z direction sum forces        
 * 
 * @property  this_node  pointer to a node so that it can access member of if 
 *                       (i.e., M, B, FX, ect.)           
 */
struct                                                   
SUM_F_IN_X                                               
{                                                        
  Node * const this_node; // can't change reference      
                          // to this_node                
  double operator()();                                   
SUM_F_IN_X( Node &nodePtr ) : this_node( &nodePtr ) {}   
};                                                       
                                                         
                                                         
struct                                                   
SUM_F_IN_Y                                               
{                                                        
  Node * const this_node; // can't change reference      
                          // to this_node                
  double operator()();                                   
SUM_F_IN_Y( Node &nodePtr ) : this_node( &nodePtr ) {}   
};                                                       
                                                         
                                                         
struct                                                   
SUM_F_IN_Z                                               
{                                                        
  Node * const this_node; // can't change reference      
                          // to this_node                
  double operator()();                                   
SUM_F_IN_Z( Node &nodePtr ) : this_node( &nodePtr ) {}   
};                                                       
//============================================================================


/** 
 * Node. Defines the forces at each element and the element connectivity. This class 
 * defines the following:
 *   - node position in X, Y and Z
 *   - external forces applied to the node FX, FY and FZ
 *   - Newton's equilibrium equation in X, Y and Z
 *
 * @property  grav                     gravitational constant [m/s^2]
 * @property  sea_density              density of seawater [kg/m^3]
 * @property  solve_X_Newton_equation  flag to identify if X direction is being solved 
 * @property  solve_Y_Newton_equation  flag to identify if Y direction is being solved 
 * @property  solve_Z_Newton_equation  flag to identify if Z direction is being solved 
 * @property  X                        X node displacement [m]
 * @property  Y                        Y node displacement [m]
 * @property  Z                        Z node displacement [m]
 * @property  M                        node mass [kg]
 * @property  B                        node displacement [m^3]
 * @property  FX                       external applied force in X [kN]
 * @property  FY                       external applied force in Y [kN]
 * @property  FZ                       external applied force in Z [kN]
 * @property  sum_FX;                  sum force in X [kN]
 * @property  sum_FY;                  sum force in Y [kN]
 * @property  sum_FZ;                  sum force in Z [kN]
 */
struct 
Node 
{
private:
  // gravitational constant and sea density, respectively
  const double grav; 
  const double sea_density;
    
  // Informs the program if the X, Y, or Z Newton equilibrium equation is solved
  // for this node.
  bool solve_X_Newton_equation; 
  bool solve_Y_Newton_equation; 
  bool solve_Z_Newton_equation; 

  // Variables that can be iterated 
  // Note: HX and HY are used for purposes of writting data to the output file    
  VarType X;  // node displacement in X global coordinates [m] 
  VarType Y;  // node displacement in Y global coordinates [m] 
  VarType Z;  // node displacement in Z global coordinates [m] 
  VarType M;  // node point mass [kg]
  VarType B;  // node displaced volumne [m^3]
  VarType FX; // node applied force in X global coordinates [kN] 
  VarType FY; // node applied force in Y global coordinates [kN]
  VarType FZ; // node applied force in Z global coordinates [kN] 
    
  // Sum forces from the contribution of all elements, applied forces
  // and mass/buoyancy modules
  double sum_FX; 
  double sum_FY; 
  double sum_FZ;         

  // Variables used for the linearization. 
  double prior_value_x;
  double prior_value_y;
  double prior_value_z;

  // return a point to this instance of a node.  
  // We use this to initialize f_x, f_y, and f_z function objects. 
  Node& self( void ) { return *this; }

public:   
  // give the two catenary equations access to the private memeber in class 
  // Element. This avoids the need for getter for all double variables in 
  // this class.
  friend class  MAP_OtherStateType_class; 
  friend class  Element;

  Node ( NodeType p , double g , double rho) : grav ( g      ) , // gravitational constant
    sea_density                                     ( rho    ) , // density of sea water
    solve_X_Newton_equation                         ( false  ) , // flag indicating the sum force equation in X is solved 
    solve_Y_Newton_equation                         ( false  ) , // flag indicating the sum force equation in Y is solved 
    solve_Z_Newton_equation                         ( false  ) , // flag indicating the sum force equation in Z is solved 
    sum_FX                                          ( 9999.9 ) , // sum forces in X direction
    sum_FY                                          ( 9999.9 ) , // sum forces in Y direction
    sum_FZ                                          ( 9999.9 ) , // sum forces in Z direction
    prior_value_x                                   ( 9999.9 ) ,
    prior_value_y                                   ( 9999.9 ) ,
    prior_value_z                                   ( 9999.9 ) , // function object to the z direction sum force equation
    type                                            ( p      ) ,
    f_x                                             ( self() ) ,
    f_y                                             ( self() ) , // function object to the y direction sum force equation 
    f_z                                             ( self() ) { }    

  ~Node ( ) { } // Kill the node
    
  // Defines the node as a Vessel, Connect or Fix type
  // @note : you cannot change the node type while the simulation is running
  const NodeType type;
    
  // f_x, f_y and f_z are instances of the force balance equation, which are used for 
  // all Connect nodes. Some equations fo f_i may be used if the node is a Vessel or
  // Fix entity.
  SUM_F_IN_X f_x;
  SUM_F_IN_Y f_y;
  SUM_F_IN_Z f_z; 

  // get data necessary for SUM_F_IN_i function objects
  double GetSumFX( ) const; 
  double GetSumFY( ) const; 
  double GetSumFZ( ) const; 
  double GetFX( ) const; 
  double GetFY( ) const; 
  double GetFZ( ) const; 
  double GetM( ) const;  
  double GetB( ) const;  
  double GetGrav( ) const; 
  double GetSeaDensity( ) const; 
        
  // Let's test if the node is fixed (an anchor).  If it isn't fixed, then let's do stuff with this node.
  template<class T>
  bool GetVarTypeFixedBool( VarType T::* ptr ) const 
    { return (this->*ptr).is_fixed; }
    
  double GetVarTypeValue( VarType Node::* ptr                ) const;
  void   SetVarTypeValue( VarType Node::* ptr , double value );
      
  NodeType GetNodeType( ) const; // get the node as a Fix, Vessel or Connect
    
  // Function used for linearization routine. The prior value for X, Y and Z is stored.
  // This is to differentiate it from x+epsilon in the finite difference process.
  void   SetEquilibriumDisplacement( );
  void   SetPriorXValue( );       
  void   SetPriorYValue( );       
  void   SetPriorZValue( );       
  double GetPriorXValue( ) const; 
  double GetPriorYValue( ) const; 
  double GetPriorZValue( ) const; 
   
  // set the 'solve_i_Newton_equation' boolean. This idicates if we are solving the
  // Newton equation to enforce static equilibrium for this node.
  void SetXNewtonEquationFlag( bool flag ); 
  void SetYNewtonEquationFlag( bool flag ); 
  void SetZNewtonEquationFlag( bool flag ); 

  // get the 'solve_i_Newton_equation' flag
  bool GetXNewtonEquationFlag( ) const;
  bool GetYNewtonEquationFlag( ) const;
  bool GetZNewtonEquationFlag( ) const;

  // Set the node either as a Vessel, Fix or Connect
  void SetVarType( const std::string &varStr , 
                   const std::string &type   , 
                   int               i       , 
                   VarType Node::*   ptr     , 
                   MAP_ErrStat_class       &err    , 
                   MAP_Message_class       &msg    );

  // initialize (or set) sum_FX, sum_FY, sum_FZ to 0.0
  void SetSumForceToZero( );

  // this sets the initial conditions for the Fix and Vessel nodes
  // This adds sum_FX, sum_FY and sum_FZ into the node
  void InitializeNode( );

  // sums the element fairlead/anchor forces to the approriate node
  void AddToSumFX( double x_force );
  void AddToSumFY( double y_force );
  void AddToSumFZ( double z_force );
};


#endif // _NODE_H
