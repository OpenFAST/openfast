/**
 * ====================================================================================================
 *                              Element.h
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


#ifndef _ELEMENT_H
#define _ELEMENT_H


#include "Node.h" /**
                   * Preprocessor Defitions in Node.h
                   * 
                   * #include "VarType.h"
                   *     #include <boost/lexical_cast.hpp>
                   *     #include <boost/algorithm/string.hpp>
                   *     #include <string>
                   *     #include <iomanip>
                   *     #include "MAP_Message.h" 
                   *     #include "MAP_ErrStat.h" 
                   */
                  
#include "CableLibrary.h" 
#include "CatenaryEquation.h"
#include <boost/smart_ptr.hpp>
#include <boost/math/special_functions/asinh.hpp>

typedef boost::shared_ptr <CableLibrary> CableLibrary_ptr;
typedef boost::shared_ptr <Node>         Node_ptr;
typedef boost::shared_ptr <Element>      Element_ptr;


/**
 * ====================================================================================================
 * Element
 *
 *
 * ====================================================================================================
 */
class Element {
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
    double psi;    // angle of roation between global X and local x axis [rad]
    double l;      // horizontal cable excursion [m]
    double h;      // vertical cable excursion [m]
    double Hx;     // horizontal fairlead force in global X direction [N]
    double Hy;     // horizontal fairlead force in global Y direction [N]
    double Vz;     // vertical anchor force [N]
    double Ha_x;   // horizontal anchor force in global X direction [N]
    double Ha_y;   // horizontal anchor force in global Y direction [N]
    double Va;     // vertical anchor force [N]
    double omega;  // cable mass per length [kg/m]
    double A;      // cross-sectional area [m^2]
    double EA;     // element stiffness [N]
    double Lb;     // length of element touching the seabed [m]
    double norm;   // distance between fairlead and anchor

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
                  EA                   ( 9999.9 ) , // element stiffness [N]
                  Lb                   ( 9999.9 ) , // length of element touching the seabed [m]
                  norm                 ( 9999.9 ) , // distance between fairlead and anchor [m]
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
    std::string getName( ) const;
    
    // Plots the cable profile in Python using the matPlotLib library
    bool elementPlotPoints( std::string &X , std::string &Y , std::string &Z , MAP_Message &Msg );

    // sets a reference to a cable property in the CableLibrary. 'line_property'
    // is a pointer
    void setLineProperty( const std::string &element_type , std::vector <CableLibrary_ptr> &T , MAP_ErrStat &Error ,  MAP_Message &Msg );
    
    // Sets a reference to the upper/lower node. 'fairlead'/'anchor' is a pointer.
    void setFairlead( const int upper , const std::vector <Node_ptr> &T );
    void setAnchor  ( const int lower , const std::vector <Node_ptr> &T );
    
    // Initialize Lu, H and V variable (name, index and value)
    void setLu( const std::string in, const int i, MAP_ErrStat &Error, MAP_Message &Msg );
    void setH ( const int index, MAP_ErrStat &Error, MAP_Message &Msg );
    void setV ( const int index, MAP_ErrStat &Error, MAP_Message &Msg );
    void setHX( ); // used solely to get access to Hx and write data to the output file
    void setHY( ); // used solely to get access to Hy and write data to the output file 

    // set H and V 'is_fixed' flags to thier default values. This is called only once at initialization. The boolean
    // values set here can be over-riden with Element::setHFlagTo() and Element::setVFlagTo(). This is called in 
    // MAP_OtherStateType_class::addElement().
    void setH_and_V_flags( );

    // set the H and V 'is_fixed' boolean manually
    void setHFlagTo( const bool T ){ this->H.is_fixed=T; }
    void setVFlagTo( const bool T ){ this->V.is_fixed=T; }

    // set element option flags to true or false
    void setPLOT_flag        ( const bool T ) { this->PLOT_flag         = T; }
    void setX_POS_flag       ( const bool T ) { this->X_POS_flag        = T; }
    void setY_POS_flag       ( const bool T ) { this->Y_POS_flag        = T; }
    void setZ_POS_flag       ( const bool T ) { this->Z_POS_flag        = T; }
    void setX_FORCE_flag     ( const bool T ) { this->X_FORCE_flag      = T; }
    void setY_FORCE_flag     ( const bool T ) { this->Y_FORCE_flag      = T; }
    void setZ_FORCE_flag     ( const bool T ) { this->Z_FORCE_flag      = T; }
    void setLINE_TENSION_flag( const bool T ) { this->LINE_TENSION_flag = T; }
    void setOMIT_CONTACT_flag( const bool T ) { this->OMIT_CONTACT_flag = T; }
    void setLAY_LENGTH_flag  ( const bool T ) { this->LAY_LENGTH_flag   = T; }
    
    // get element option flags to true or false
    bool getPLOT_flag        ( ) const { return this->PLOT_flag;         }
    bool getX_POS_flag       ( ) const { return this->X_POS_flag;        }
    bool getY_POS_flag       ( ) const { return this->Y_POS_flag;        }
    bool getZ_POS_flag       ( ) const { return this->Z_POS_flag;        }
    bool getX_FORCE_flag     ( ) const { return this->X_FORCE_flag;      }
    bool getY_FORCE_flag     ( ) const { return this->Y_FORCE_flag;      }
    bool getZ_FORCE_flag     ( ) const { return this->Z_FORCE_flag;      }
    bool getLINE_TENSION_flag( ) const { return this->LINE_TENSION_flag; }
    bool getOMIT_CONTACT_flag( ) const { return this->OMIT_CONTACT_flag; }
    bool getLAY_LENGTH_flag  ( ) const { return this->LAY_LENGTH_flag;   }

    // returns the VarType::value parameter
    double getLu( ) const ;
    double getH ( ) const ; 
    double getV ( ) const ;
 
    // print line tension to map output file
    std::string getLineTension      (             );
    std::string getLineTensionHeader( const int i );
    
    // returns the VarType::is_fixed parameter
    bool getLuFlag( ) const;
    bool getHFlag ( ) const;
    bool getVFlag ( ) const;
     
    // returns the Node::X.value, Node::Y.value, Node::Z.value for the 
    // fairlead/anchor
    double getAnchorPosition  ( VarType Node::* ptr ) const;
    double getFairleadPosition( VarType Node::* ptr ) const;

    // return the xf, xa, yf, ya value for each node. This is used in A_DERIVS in UserData.h
    double getXf( ) { return this->getFairleadPosition( &Node::X ); }
    double getXa( ) { return this->getAnchorPosition  ( &Node::X ); }
    double getYf( ) { return this->getFairleadPosition( &Node::Y ); }
    double getYa( ) { return this->getAnchorPosition  ( &Node::Y ); }

    // returns the Node::X.is_fixed, Node::Y.is_fixed, Node::Z.is_fixed 
    // for the fairlead/anchor    
    bool getAnchorFlag  ( VarType Node::* ptr ) const;
    bool getFairleadFlag( VarType Node::* ptr ) const;
    
    // set element variables
    void initializeElement( const double &g, const double &rho , MAP_ErrStat &err , MAP_Message &Msg ); // throws
    void updateElement    ( MAP_ErrStat &err , MAP_Message &Msg );

    // the following four member functions are called at element
    // initialization
    void   setPsi                 (                                     ); // throws
    double getPsi                 (                                     ) const { return this->psi; }
    void   initializeCableProperty( const double &g , const double &rho ); // throws
    void   set_l_and_h            (                                     );
    void   initializeHAndV        ( const double &g , const double &rho );
    
    // set the variables 'sum_FX', 'sum_FY' and 'sum_FZ' for the 'fairlead' and 
    // 'anchor' nodes to zero. This method calls Node::setSumForceToZero()
    void cleanNodes( );

    // We add the negative of the element anchor force to the node 
    // (Node::sum_FX, Node::sum_FY, Node::sum_FY)
    void addForceToFairleadNode( ); 
    void addForceToAnchorNode  ( );    
    
    // get value of VarType::H.value, VarType::V.value
    //double getElementVar   ( VarType Element::* ptr ) { return (this->*ptr).value;  }
    double getOmega        (                 ) const { return this->omega;                   }
    double getHeight       (                 ) const { return this->h;                       }
    double getLength       (                 ) const { return this->l;                       }
    double getArea         (                 ) const { return this->A;                       }
    double getEA           (                 ) const { return this->line_property->EA.value; }
    double getCB           (                 ) const { return this->line_property->CB.value; }
    bool   getRestingSeabed(                 ) const { return is_resting_on_seabed;          }
    void   setRestingSeabed( const bool flag )       { this->is_resting_on_seabed = flag;    }

    bool compareNodeAddressWithFairlead ( const Node &ref ) const;
    bool compareNodeAddressWithAnchor   ( const Node &ref ) const;
};


#endif // _ELEMENT_H
