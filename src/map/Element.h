/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   Element.h
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


#ifndef _ELEMENT_H
#define _ELEMENT_H


#include "Node.h" // Preprocessor Defitions in Node.h
                  //   #include "VarType.h"
                  //       #include <boost/lexical_cast.hpp>
                  //       #include <boost/algorithm/string.hpp>
                  //       #include <string>
                  //       #include <iomanip>
                  //       #include "MAP_Message_class.h" 
                  //       #include "MAP_ErrStat_class.h" 
                  
#include "CableLibrary.h" 
#include "CatenaryEquation.h"
#include <boost/smart_ptr.hpp>
#include <boost/math/special_functions/asinh.hpp>

typedef boost::shared_ptr <CableLibrary> CableLibrary_ptr;
typedef boost::shared_ptr <Node>         Node_ptr;
typedef boost::shared_ptr <Element>      Element_ptr;


// ====================================================================================================
// Element
//
// ====================================================================================================
class 
Element 
{
private:

  // line properties
  CableLibrary *line_property;
    
  // informs if the cable is resting on the seabed so the plot member funtion in 
  // MAP_OtherStateType_class is aware of the equations to used. 
  bool is_resting_on_seabed;
    
  // anchor and fairlead node, respectively
  Node *anchor;
  Node *fairlead;

  MAP_ERROR_CODE error_code;

  // line/element variables
  double psi;           // angle of roation between global X and local x axis [rad]
  double l;             // horizontal cable excursion [m]
  double h;             // vertical cable excursion [m]
  double Hx;            // horizontal fairlead force in global X direction [kN]
  double Hy;            // horizontal fairlead force in global Y direction [kN]
  double Vz;            // vertical fairlead force [kN]
  double Ha_x;          // horizontal anchor force in global X direction [kN]
  double Ha_y;          // horizontal anchor force in global Y direction [kN]
  double Va;            // vertical anchor force [kN]
  double omega;         // cable weight per length in seawater [N/m]
  double A;             // cross-sectional area [m^2]
  double EA;            // element stiffness [kN]
  double Lb;            // length of element touching the seabed [m]
  double norm;          // distance between fairlead and anchor
  double ghost_value_x; // This is used for the case of a vertical cable. ='s original x fairlead value
  double ghost_value_y; // This is used for the case of a vertical cable. ='s original y fairlead value
  bool   is_vertical;   // is the cable vertical? True = yes

  // Element option flags from the MAP input file
  bool PLOT_flag; 
  bool X_POS_flag; 
  bool Y_POS_flag; 
  bool Z_POS_flag;
  bool X_FORCE_flag; 
  bool Y_FORCE_flag; 
  bool Z_FORCE_flag;
  bool LINE_TENSION_flag; 
  bool OMIT_CONTACT_flag; 
  bool LAY_LENGTH_flag; 
  bool FAIR_TENSION_flag;
  bool ANCH_TENSION_flag;

  // A vector array of line tension along the length of the cable.
  // 
  // @todo : current, the number of points we evaluate the tension is fixed to 10.
  //         This should be changed in a future version of MAP so the number of 
  //         points evaluated is selected at run time. 
  // 
  std::vector <double> line_tension;
    
  // Variables that can be iterated
  // Note: HX and HY are used for purposes of writting data to the output file
  VarType Lu; // unstretched cable length
  VarType H;  // Horizontal fairlead force in the local cable elemenet frame
  VarType V;  // Vertical fairlead force
  VarType HX; 
  VarType HY;      

  // return a point to this instance of a node.  
  // We use this to initialize f_v and f_h function objects. 
  Element &self( void ) { return *this; }
    
public:    
  // give the two catenary equations access to the private memeber in class 
  // Element. This avoids the need for getter for all double variables in 
  // this class.
  friend class MAP_OtherStateType_class; 

  Element ( ) : line_property        ( NULL   ) , // pointer to the CableLibrary class 
    is_resting_on_seabed ( false  ) , // flag : identifying if the cable is on the sea floor
    anchor               ( NULL   ) , // pointer to anchor node
    fairlead             ( NULL   ) , // pointer to fairlead node
    error_code           (MAP_SAFE) ,
    psi                  ( 9999.9 ) , // relative angle between the element x and global X axis
    l                    ( 9999.9 ) , // cable span in xy plane
    h                    ( 9999.9 ) , // cable height; equal to (fairlead.z)-(anchor.z)
    Hx                   ( 9999.9 ) , // x fairlead force in the global frame
    Hy                   ( 9999.9 ) , // y fairlead force in the global frame
    Vz                   ( 9999.9 ) , // z fairlead force in both the global and element frame
    Ha_x                 ( 9999.9 ) , // x anchor force in the global frame
    Ha_y                 ( 9999.9 ) , // y anchor force in the global frame
    Va                   ( 9999.9 ) , // z anchor force in both the global and element frame
    omega                ( 9999.9 ) , // cable weight per unit length in seawater
    A                    ( 9999.9 ) , // cable cross-sectional area [m^2] 
    EA                   ( 9999.9 ) , // element stiffness [kN]
    Lb                   ( 9999.9 ) , // length of element touching the seabed [m]
    norm                 ( 9999.9 ) , // distance between fairlead and anchor [m]
    ghost_value_x        ( 0.0    ) , // ='s original x fairlead value
    ghost_value_y        ( 0.0    ) , // ='s original y fairlead value
    is_vertical          ( false  ) , // is the cable vertical? True = yes
    PLOT_flag            ( false  ) , // flag : plot this element only 
    X_POS_flag           ( false  ) , // flag : write global x position data to text file
    Y_POS_flag           ( false  ) , // flag : write global y position data to text file
    Z_POS_flag           ( false  ) , // flag : write global z position data to text file
    X_FORCE_flag         ( false  ) , // flag : write global x force data to text file
    Y_FORCE_flag         ( false  ) , // flag : write global y force data to text file
    Z_FORCE_flag         ( false  ) , // flag : write global z force data to text file
    LINE_TENSION_flag    ( false  ) , // flag : write line tension magnitude to text file
    OMIT_CONTACT_flag    ( false  ) , // flag : neglect contact with seabed
    LAY_LENGTH_flag      ( false  ) , // flag : length of line laying on the sea floor
    FAIR_TENSION_flag    ( false  ) , // flag : fairlead node tension
    ANCH_TENSION_flag    ( false  ) , // flag : anchor node tension
    f_h                  ( self() ) , // function object to the horizontal catenary equation
    f_v                  ( self() ) , 
    dXdH                 ( self() ) , 
    dXdV                 ( self() ) , 
    dZdH                 ( self() ) , 
    dZdV                 ( self() ) {} // function object to the vertical catenary equation

  ~Element ( ) { } // destructor

  // function objects to the two catenary equations
  HORIZONTAL_CATENARY_EQ f_h;
  VERTICAL_CATENARY_EQ   f_v;
  DXDH                   dXdH;
  DXDV                   dXdV;
  DZDH                   dZdH;
  DZDV                   dZdV;    

  // returns the type of cable which defines this element, 
  // i.e., nylon, steel, ect.
  std::string GetElementName( ) const;
    
  // Plots the cable profile in Python using the matPlotLib library
  bool GetElementPlotArray( std::string &X, std::string &Y, std::string &Z, MAP_ErrStat_class &err, MAP_Message_class &Msg );
  void GetElementPlotArray( std::vector<double> &X, std::vector<double> &Y, std::vector<double> &Z, MAP_ErrStat_class &err, MAP_Message_class &msg );

  // sets a reference to a cable property in the CableLibrary. 'line_property'
  // is a pointer
  void SetLineProperty( const std::string              &element_type , 
                        std::vector <CableLibrary_ptr> &library      , 
                        MAP_ErrStat_class                    &err          ,  
                        MAP_Message_class                    &msg          );
    
  // Sets a reference to the upper/lower node. 'fairlead'/'anchor' is a pointer.
  void SetFairlead( const int upper , const std::vector <Node_ptr> &node_library );
  void SetAnchor( const int lower , const std::vector <Node_ptr> &node_library );
    
  // Initialize Lu, H and V variable (name, index and value)
  void SetV( const int index, MAP_ErrStat_class &err, MAP_Message_class &msg );
  void SetH( const int index, MAP_ErrStat_class &err, MAP_Message_class &msg );
  void SetLu( const std::string  varStr, const int i, MAP_ErrStat_class &err, MAP_Message_class &msg );

  void SetHX( ); // used solely to get access to Hx and write data to the output file
  void SetHY( ); // used solely to get access to Hy and write data to the output file 

  // set H and V 'is_fixed' flags to thier default values. This is called only once at initialization. The boolean
  // values set here can be over-riden with Element::SetHFlagTo() and Element::SetVFlagTo(). This is called in 
  // MAP_OtherStateType_class::addElement().
  void SetH_and_V_flags( );
  void CheckIfCableIsVertical();
  void ResetGhostProperties();

  // set the H and V 'is_fixed' boolean manually
  void SetHFlagTo( const bool flag );
  void SetVFlagTo( const bool flag );

  // set element option flags to true or false
  void SetPlotFlag( const bool flag ); 
  void SetXPosFlag( const bool flag );
  void SetYPosFlag( const bool flag );
  void SetZPosFlag( const bool flag );
  void SetXForceFlag( const bool flag );
  void SetyForceFlag( const bool flag );
  void SetZForceFlag( const bool flag );
  void SetLineTensionFlag( const bool flag );
  void SetOmitContactFlag( const bool flag );
  void SetLayLengthFlag  ( const bool flag );
  void SetFairTensionFlag( const bool flag );
  void SetAnchTensionFlag( const bool flag );
    
  // get element option flags to true or false
  bool GetPlotFlag( ) const; 
  bool GetXPosFlag( ) const;
  bool GetYPosFlag( ) const;
  bool GetZPosFlag( ) const;
  bool GetXForceFlag( ) const;
  bool GetYForceFlag( ) const;
  bool GetZForceFlag( ) const;
  bool GetLineTensionFlag( ) const;
  bool GetOmitContactFlag( ) const;
  bool GetLayLengthFlag( ) const;

  // returns the VarType::value parameter
  double GetLu( ) const ;
  double GetH( ) const ; 
  double GetV( ) const ;
 
  // print line tension to map output file
  std::string GetLineTensionString( );
  std::string GetLineTensionStringHeader( const int i );
  std::string GetLineTensionStringHeaderForFast( const int i, const int j );
  std::string GetLayLengthString( );  
  std::string GetLayLengthStringHeader( const int i );  
  std::string GetFairleadTensionStringHeader( const int i );
  std::string GetAnchorTensionStringHeader( const int i );

  double GetFairleadTensionMagnitude( ) const;
  double GetAnchorTensionMagnitude( ) const;
  std::string GetFairleadTensionMagnitudeString( ) const;
  std::string GetAnchorTensionMagnitudeString( ) const;

  // returns the VarType::is_fixed parameter
  bool  GetLuFlag( ) const;
  bool  GetHFlag( ) const;
  bool  GetVFlag( ) const;
  float GetLayLengthValue() const;
  void  GetLineTensionValues( float * passArr );

  // returns the Node::X.value, Node::Y.value, Node::Z.value for the 
  // fairlead/anchor
  double GetAnchorPosition  ( VarType Node::* ptr ) const;
  double GetFairleadPostion( VarType Node::* ptr ) const;

  // return the xf, xa, yf, ya value for each node. This is used in A_DERIVS in UserData.h
  double GetXf( ) const;
  double GetXa( ) const;
  double GetYf( ) const;
  double GetYa( ) const;

  // returns the Node::X.is_fixed, Node::Y.is_fixed, Node::Z.is_fixed 
  // for the fairlead/anchor    
  bool GetAnchorFlag( VarType Node::* ptr ) const;
  bool GetFairleadFlag( VarType Node::* ptr ) const;
    
  // set element variables
  void InitializeElement( const double &g, const double &rho_sea, MAP_ErrStat_class  &err, MAP_Message_class  &msg ); // throws MAP_ERROR_56 and/or MAP_WARNING_1

  void UpdateElement( MAP_ErrStat_class &err , MAP_Message_class &msg );

  // the following four member functions are called at element initialization
  void SetCableDistance( );
  void CheckMaximumLineLength( ); // throws
  void SetPsi( ); // throws
  double GetPsi( ) const;
  void InitializeCableProperty( const double &g , const double &rho_sea ); // throws MAP_WARNING/MAP_EROR_56
  void InitializeHAndV( const double &g , const double &rho_sea );
    
  // set the variables 'sum_FX', 'sum_FY' and 'sum_FZ' for the 'fairlead' and 
  // 'anchor' nodes to zero. This method calls Node::SetSumForceToZero()
  void ResetNodes( );

  // We add the negative of the element anchor force to the node 
  // (Node::sum_FX, Node::sum_FY, Node::sum_FY)
  void AddForceToFairleadNode( ); 
  void AddForceToAnchorNode( );    
    
  // get value of VarType::H.value, VarType::V.value
  //double getElementVar   ( VarType Element::* ptr ) { return (this->*ptr).value;  }
  double GetOmega( ) const; 
  double GetHeight( ) const; 
  double GetLength( ) const; 
  double GetArea( ) const; 
  double GetEA( ) const; 
  double GetCB( ) const; 
  bool GetRestingSeabedFlag( ) const; 
  void SetRestingSeabedFlag( const bool flag );

  bool CompareNodeAddressWithFairlead( const Node &ref ) const;
  bool CompareNodeAddressWithAnchor( const Node &ref ) const;
};


#endif // _ELEMENT_H
