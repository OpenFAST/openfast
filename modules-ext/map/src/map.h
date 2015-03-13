/****************************************************************
 *   Copyright (C) 2014 mdm                                     *
 *   map[dot]plus[dot]plus[dot]help[at]gmail                    *
 *                                                              *
 * Licensed to the Apache Software Foundation (ASF) under one   *
 * or more contributor license agreements.  See the NOTICE file *
 * distributed with this work for additional information        *
 * regarding copyright ownership.  The ASF licenses this file   *
 * to you under the Apache License, Version 2.0 (the            *
 * "License"); you may not use this file except in compliance   *
 * with the License.  You may obtain a copy of the License at   *
 *                                                              *
 *   http://www.apache.org/licenses/LICENSE-2.0                 *
 *                                                              *
 * Unless required by applicable law or agreed to in writing,   *
 * software distributed under the License is distributed on an  *
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       *
 * KIND, either express or implied.  See the License for the    *
 * specific language governing permissions and limitations      *      
 * under the License.                                           *  
 ****************************************************************/


#ifndef _MAP_H
#define _MAP_H


#include "mapsys.h"

#include "simclist/simclist.h"

#include "bstring/bstrlib.h"

#include "cminpack/cminpack.h"
#include "cminpack/cminpackP.h"
#include "cminpack/minpack.h"

#include "MAP_Types.h"
#include "maperror.h"
#include "lmroutines.hpp"

#ifdef WITH_LAPACK
#  include "lapack/lapacke.h"
#endif // WITH_LAPACK

/**
 * @brief Associates the node with a particular type. Fix nodes are anchor points
 *        (ussually) and cannot move with time. Connect nodes are intermediaries
 *        connecting two lines. When a connect node is defined, the force-balance
 *        equation is minimized to obtain the cable profile X, Y, and Z node position
 *        at connect nodes are associated with constraint types. Vessel nodes attach
 *        to the platform. X, Y, and Z, node positions are associated with input 
 *        types, and the corresponding fx, fy, and fz forces are the outputs.  
 *        Moments are calculated by the calling program. 
 */
typedef enum NodeType_enum {
  NONE,     /**< */
  FIX,      /**< */
  CONNECT,  /**< */
  VESSEL    /**< */
} NodeType;


/**
 * @brief Finite different routine when solving the outer-loop iterations. 
 *        This is also used when the linearization routine is called. 
 */
typedef enum FdType_enum {
  BACKWARD_DIFFERENCE, /**< */
  CENTRAL_DIFFERENCE,  /**< */
  FORWARD_DIFFERENCE   /**< */
} FdType;


/**
 * @struct SolveType
 * @brief Set internally. If a connect node is present in the MAP input
 *        file, then the solve typeis partitioned. Partitioned describes the
 *        process of breaking the problem into two blocks: an inner-loop solve
 *        (where the catenary equations are minimized) and an outer-loop solve
 *        (where the force-balance equations are minimized). Likewise, when no
 *        connect nodes are present, the solver type is monolithic and only the
 *        non-linear catenary equations need to be solved. 
 */
typedef enum SolveType_enum {  
  MONOLITHIC,   /**< for MSQS, lines only (no connect nodes) */
  PARTITIONED,  /**< for MSQS system with connect nodes */
  LUMPED_MASS   /**< for lumped-mass model */
} SolveType;


/**
 * @brief Finite-difference structure. This is used locally withint FD routines to conveniently store the
 *        force and moment when vessel is displaced by epsilon. 
 */
struct Fd_t {
  double* fx;
  double* fy;
  double* fz;
  double* mx;
  double* my;
  double* mz;
}; typedef struct Fd_t Fd;


/**
 * @brief Fundamental MAP type. Unless it is used locally, every variable should be darclared as a VarType
 *        (or VarTypePtr, described below). The VarType provides convenience for printing variable 
 *        information, such as the units, it's name (such as 'FX[1]' for X-direction node 1 force to the 
 *        output file), and counting references.
 */
struct VarType_t {
  bstring units;    /**< units for printing information to a summary file or output buffer */
  bstring name;     /**< name of the variable. This is used for identifying it in the output buffer */
  double value;     /**< the value */
  bool is_fixed;    /**< if is_fixed = true, then we are not solving for this variable */
  bool user_initial_guess; /**< if user_initial_guess = true, the user has supplied an initial guess */
  int ref_counter;  /**< for ensuring the variable is assigned to one of: input, param, or constraint */
  int id;           /**< node or line this value is attached to */
}; typedef struct VarType_t VarType;


/**
 * @brief Serves the same function as VarType, but treats value as a pointer. This preserves the FAST integration
 *        to native Fortran derivived types. Instead of the variable residing in C, the value parameter points
 *        to a variable allocated in Fortran. This feature is also preserved with Python binding. 
 */
struct VarTypePtr_t {
  bstring units;    /**< units for printing information to a summary file or output buffer */
  bstring name;     /**< name of the variable. This is used for identifying it in the output buffer */
  double* value;   /**< the value */
  bool is_fixed;    /**< If is_fixed = true, then we are not solving for this variable */
  int ref_counter;  /**< For ensuring the variable is assigned to one of: input, param, or constraint */
  int id;           /**< node or line this value is attached to */
}; typedef struct VarTypePtr_t VarTypePtr;


struct Vector_t {
  double x;
  double y;
  double z;
}; typedef struct Vector_t Vector;
  

struct Point_t {
  VarType x;
  VarType y;
  VarType z;  
}; typedef struct Point_t Point;


struct PointPtr_t {
  VarTypePtr x;
  VarTypePtr y;
  VarTypePtr z;  
}; typedef struct PointPtr_t PointPtr;


struct EulerAngle_t {
  VarType phi;
  VarType the;
  VarType psi;  
}; typedef struct EulerAngle_t EulerAngle;


struct Force_t {
  VarType fx;
  VarType fy;
  VarType fz;  
}; typedef struct Force_t Force;


struct ForcePtr_t {
  VarTypePtr fx;
  VarTypePtr fy;
  VarTypePtr fz;  
}; typedef struct ForcePtr_t ForcePtr;


struct ReferencePoint_t {
  VarTypePtr* x;
  VarTypePtr* y;
  VarTypePtr* z;
}; typedef struct ReferencePoint_t ReferencePoint;


/**
 * @brief Central point where all 'VESSEL' nodes can be displaced. Instead of displacing all nodes individually, the vessel can be displaced,
 *        then helper functions can be called to displacement the nodes. The vessel only reference 'input' nodes.
 */
struct Vessel_t {
  EulerAngle orientation; /**< Vessel orientation [deg]*/
  double* xi;             /**< initial node connection point in body frame.This is equal to uType->x at initialization [m] */
  double* yi;             /**< initial node connection point in body frame This is equal to uType->y at initialization [m] */
  double* zi;             /**< initial node connection point in body frame This is equal to uType->z at initialization [m] */
  Point displacement;     /**< User-specified vessel displacement. This is the [m]*/
  Point ref_origin;       /**< Center of rotation origin. The moments are taken about this is point. The reference point is with respect to the FAST reference origin (equal to the
                           * SWL at zero vessel dispalcements) [m]*/
  Force line_sum_force;   /**< Sum force of all nodes connecting to the vessel [N] */
}; typedef struct Vessel_t Vessel;



/**
 * @brief References a list of VarType's (in the case of out_list) and VarTypePtr's
 *        (in the case of out_list_ptr) for output file data streaming. The idea is
 *        to loop through the list and print the 'value' and 'name[id]' members of
 *        VarType and VarTypePtr.
 * @see   
 * @todo  This really should be made polymorphic so only one variable is needed. 
 */
struct OutputList_t{
  list_t out_list;     /**< Outputs associated with VarType */
  list_t out_list_ptr; /**< Outputs associated with VarTypePtr */
}; typedef struct OutputList_t OutputList;


typedef struct {
  bool gx_pos_flag;
  bool gy_pos_flag;
  bool gz_pos_flag;
  bool gx_anchor_pos_flag;
  bool gy_anchor_pos_flag;
  bool gz_anchor_pos_flag;
  bool gx_force_flag;
  bool gy_force_flag;
  bool gz_force_flag;
  bool H_flag; 
  bool V_flag;
  bool V_anchor_flag;
  bool H_anchor_flag;
  bool fairlead_tension_flag;
  bool anchor_tension_flag;
  bool horizontal_excursion_flag;
  bool vertical_excursion_flag;
  bool azimuth_flag;
  bool altitude_flag;
  bool altitude_anchor_flag;
  bool line_tension_flag;
  bool omit_contact;
  bool lay_length_flag;
  bool damage_time_flag;
  bool diagnostics_flag;
  bool linear_spring;                  /**< treat the elastic catenary as a uncompressible linear springs when true */         
} LineOptions; 


/**
 * @brief Defines cable properties for a line. These values are fixed with time and connot change.
 */
struct CableLibrary_t {
  double diam;         /**< Cable diameter, [m] */
  double mass_density; /**< Cable density in air [kg/m] */
  double EA;           /**< Line stiffness [N] */
  double omega;        /**< cable weight per length in seawater [N/m] */
  double a;            /**< cross-sectional area [m^2] */
  double cb;           /**< Cable/seabed friction coefficient [non-dimensional] */
  double cd_i;         /**< Internal (structural) damping coefficient [non-dimensional] */
  double ca;           /**< Added mass coefficient [non-dimensional] */
  double cd_n;         /**< Quadtradice drag coefficient in the cable cross-flow direction [non-dimensional] */
  double cd_t;         /**< Tangential drag oefficient [non-dimensional] */
  bstring label;        /**< Provides the string a recognizable name (such as 'nylon' or 'steel') */
}; typedef struct CableLibrary_t CableLibrary;


struct Node_t {
  PointPtr acceleration;        /**< Node accelration; integrated quantity [m/s^2]; used for LM model. Associated with continuous type */
  PointPtr velocity;            /**< Node velocity [m/s]; used for LM model. Associated with continuous type */
  PointPtr position_ptr;        /**< this is a Ptr because it points to a fortran type */
  NodeType type;
  ForcePtr sum_force_ptr;       /**< this is a Ptr because it points to a fortran type */ 
  VarType M_applied;
  VarType B_applied;
  Force external_force;    
}; typedef struct Node_t Node;


// struct LMAttributes {
//   void* lm_container; /**< container struct in lmroutines.hpp */
//   double* FlineS;     /**< retains last solution for when this is called with dT = 0 */   // 
//   double** rFairtS;   /**< fairlead locations ON TURBINE */                               // <--------- will be mapped to a Node_t list 
//   double** rFairRel;  /**< fairlead locations relative to platform center */              // <--------- will be mapped to a Node_t list 
//   double** rFairi;    /**< inertial Fairlead Locations  */                                // <--------- will be mapped to a Node_t list    
//   double** rdFairi;   /**< inertial Fairlead Velocities */                                // <--------- will be mapped to a Node_t list  
// 
//   // static vectors to hold line and connection objects!
//   // vector< LineProps > LinePropList; // to hold line library types   <--------- Moved to the Line_t struct
//   // vector< Line > LineList;          //  global and persistent?      <--------- Moved to a container struct in lmroutines.hpp to isolate c++ and c
//   // vector< Connection > ConnectList;                                 <--------- Moved to a container struct in lmroutines.hpp to isolate c++ and c
//   int nConnects; 
//   int nLines;
//   
//   // state vector and stuff
//   double* states;     /**< pointer to array comprising global state vector */
//   double* newstates;
//   int Nx;             /**< size of state vector array */
// 
//   // more state vector things for rk4 integration 
//   double* f0;
//   double* f1;
//   double* f2;
//   double* f3;
//   double* xt; 
// 
//   int* lineStateIs; /** vector of line starting indices in "states" array */ // vector< int > LineStateIs;  
//   int closed; // initialize to 0
//   double dt;  // FAST time step, @rm
//   double dts; /**< mooring line time step */
// }; typedef struct LMAttributes_t LMAttributes;


struct Line_t {
  CableLibrary* line_property; /**< line properties */
  // LMAttributes* lm_attributes; /**< Preserves information of the LM model from previous time steps */
  LineOptions options;         /**< run-time options flag */
  VarTypePtr H;                /**< Horizontal fairlead force in the local cable elemenet frame */
  VarTypePtr V;                /**< Vertical fairlead force in the local cable elemenet frame */  
  VarType Lu;                  /**< unstretched cable length [m] */
  bstring label;               /**< reference a pre-defined property in the line dictionary */
  double* line_tension;        /**< array of line tension along 's' [N] */ 
  double psi;                  /**< angle of roation between global X and local x axis [deg] */
  double alpha;                /**< angle of inclication [deg] */
  double alpha_at_anchor;      /**< angle of inclication at anchor [deg] */
  double l;                    /**< horizontal cable excursion [m] */
  double Lb;                   /**< length of line touching the seabed [m] */
  double h;                    /**< vertical cable excursion [m] */
  double H_at_anchor;          /**< Horizontal anchor force in the local cable elemenet frame */
  double V_at_anchor;          /**< Vertical anchor force in the local cable elemenet frame */
  double T;                    /**< Tension magnitude [N] */
  double T_at_anchor;          /**< Tension magnitude at anchor [N] */
  double damage_time;          /**< time to damage this line and return zero force to the glue code */
  double residual_norm;       
  double fx_fairlead;
  double fy_fairlead;
  double fz_fairlead;
  double fx_anchor;
  double fy_anchor;
  double fz_anchor;
  list_t elements;             /**< LM model elements */
  Node* anchor;                /**< Anchor node */
  Node* fairlead;              /**< Fairlead node */
  int segment_size;
  int diagnostic_type;         /**< none=0, first iteration only=2, all iterations otherwise */
  int evals;                   /**< number of function evaluations */ 
  int njac_evals;              /**< number of function evaluations */      
  int converge_reason;         /*   - info=0 : improper input parameters.
                                *   - info=1 : both actual and predicted relative reductions in the sum of squares are at most ftol.
                                *   - info=2 : relative error between two consecutive iterates is at most xtol.
                                *   - info=3 : conditions for info = 1 and info = 2 both hold.
                                *   - info=4 : the cosine of the angle between fvec and any column of the Jacobian is at most gtol in absolute value.
                                *   - info=5 : number of calls to fcn has reached or exceeded maxfev.
                                *   - info=6 : ftol is too small. No further reduction in the sum of squares is possible.
                                *   - info=7 : xtol is too small. No further improvement in the approximate solution x is possible.
                                *   - info=8 : gtol is too small. fvec is orthogonal to the columns of the Jacobian to machine precision.
                                */ 
}; typedef struct Line_t Line;


struct Element_t {
  double l;  /**< \left \| \mathbf{r}_{1}-\mathbf{r}_{2} \right \| */
  Node* r1;  /**< upper node */
  Node* r2;  /**< lower node */
}; typedef struct Element_t Element;


struct InnerSolveAttributes_t {
  double** node_jac;
  double f_tol; 
  double g_tol; 
  double x_tol; 
  double x[2];               /**< array of variables the length of n */
  double fvec[2];            /**< function evaluations (residual) */
  double fjac[4];            /**< jacobian. This is a little convoluted because the jacobian is not an array */
  double wa1[2];             /**< work array of length n */
  double wa2[2];             /**< work array of length n */
  double wa3[2];             /**< work array of length n */
  double wa4[2];             /**< work array of length m */  
  double diag[2];            
  double qtf[2];             
  double factor;             
  int max_its;   
  int ldfjac;                 /**< number of columns in fjac */
  int mode;             
  int nprint;           
  int info;             
  int ipvt[2];
  int m;                      /**< number of functions */ 
  int n;                      /**< number of variables */
}; typedef struct InnerSolveAttributes_t InnerSolveAttributes;



struct OuterSolveAttributes_t {
  FdType fd;
  double** AV;        /**< for the Krylov accelerator */
  double** V;         /**< for the Krylov accelerator */
  double* av;         /**< for the Krylov accelerator, rown-major storage for AV */
  double** jac;
  double** l;         /**< lower triangle matrix in LU */
  double** u;         /**< upper triangle matrix in LU */
  double* b;          /**< this is the force vector used in x += ([J]^-1)*b */
  double* C;          /**< for the Krylov accelerator */
  double* w;          /**< for the Krylov accelerator */
  double* q;          /**< for the Krylov accelerator */
  double* x;
  double* y;
  double ds;
  double d;
  double epsilon;
  double norm_error;
  double tol;
  double coef;
  bool pg;
  bool krylov_accelerator;
  bool powell;
  int max_krylov_its;
  int max_its;
  int iteration_count;
}; typedef struct OuterSolveAttributes_t OuterSolveAttributes;


struct DomainOptions_t {
  double* repeat_angle;
  double integration_dt;  /**< Integration time step [sec]. LM model specific */
  double kb_lm;           /**< Seabed stiffeness coefficient [N/m]. LM model specific */
  double cb_lm;           /**< Seabed damping parameter [N-s/m]. LM model specific */
  bool wave_kinematics;   /**< Enable wave kinematics o calculated relative flui velcity. LM model specific */
  bool lm_model;          /**< use the lumped-mass model when true */         
  int repeat_angle_size;
}; typedef struct DomainOptions_t DomainOptions;


struct Domain_t {
  InnerSolveAttributes inner_loop; /**< Inner-loop (line level) solver options. Default settings in {@link initialize_solver_data_to_null} */
  OuterSolveAttributes outer_loop; /**< Outer-loop (node level) solver options. Default settings in {@link initialize_model_options_to_defaults} */  
  DomainOptions model_options;     /**< Contains global model options. Default setting in {@link initialize_model_options_to_defaults} */
  OutputList* y_list;              /**< Output stream. Set to null at initialization */
  SolveType MAP_SOLVE_TYPE;        /**< Identifies the solver type: single line, partitioned (multisegmented), and lumped-mass/FEA. Initialized in {@link initialize_domain_to_null}
                                    * */
  Vessel vessel;                   /**< Vessel for the mooring instance. Initialized in {@link initialize_vessel_to_null}. Associated VarType's are set in {@link set_vessel} */
  list_t library;                  /**< Cable library link list; stores cable properties, e.g., @see CableLibrary_t */
  list_t line;                     /**< Line link list */
  list_t node;                     /**< Node link list */
  list_t u_update_list;            /**< List to update the references in VarType-associated u_type's in UpdateStates. Used when coupled to FAST */
  void* HEAD_U_TYPE;               /**< Checks if the reference to MAP_InputType_t changes */
}; typedef struct Domain_t Domain;


/**
 * @details MAP options from parsed input file. Note that MAP does not readon the input file. This is done by the calling program.
 *          The calling program simply sets library_input_string, node_input_string, line_input_string, and solver_options_string.
 *          MAP then parses this string and expands them if necessary depending on the '{@link DomainOptions_t}' repeat_angle flag.
 */
struct InitializationData_t {
  struct bstrList* library_input_string;          /**< library property string from input file. MAP does not read contents from input string; must be done by calling program */
  struct bstrList* node_input_string;             /**< raw (non-expanded) node input string. MAP does not read contents from input string; must be done by calling program */
  struct bstrList* line_input_string;          /**< raw (non-expanded) line input string(MAP does not read contents from input string; must be done by calling program */
  struct bstrList* solver_options_string;         /**< model poptions input string */
  struct bstrList* expanded_node_input_string;    /**< full node input string duplicating information in node_input_string when the 'repeat' flag is used */
  struct bstrList* expanded_line_input_string; /**< full line input string duplicating information in nodeLineString when the 'repeat' flag is used */
  bstring summary_file_name;                      /**< summary file name. Can be set through {@link map_set_summary_file_name()} */
}; typedef struct InitializationData_t InitializationData;


#endif /* _MAP_H */
