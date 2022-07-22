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


#ifndef _MAPINIT_H
#define _MAPINIT_H


#include "map.h"


MAP_ERROR_CODE initialize_fortran_types(MAP_InputType_t* u_type, MAP_ParameterType_t* p_type, MAP_ContinuousStateType_t* x_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, MAP_InitOutputType_t* initout_type);
MAP_ERROR_CODE allocate_outlist(Domain* data, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Returns the size of CableLibrary structure for link list creation
 * @param   el, opaque object used in simclist
 * @see     map_init()
 * @return  Size of CableLibrary structure
 */
size_t cable_library_meter(const void* el);


/**
 * @brief   Returns the size of Node structure for link list creation
 * @param   el, opaque object used in simclist
 * @see     map_init()
 * @return  Size of Node structure
 */
size_t node_meter(const void* el);


/**
 * @brief   Returns the size of Line structure for link list creation
 * @param   el, opaque object used in simclist
 * @see     map_init()
 * @return  Size of Line structure
 */
size_t cable_line_meter(const void* el);


/**
 * @brief   Returns the size of ReferencePoint (inputs) structure for link list creation
 * @param   el, opaque object used in simclist
 * @see     map_init()
 * @return  Size of ReferencePoint (inputs) structure
 */
size_t u_list_meter(const void *el); 


/**
 * @brief   Initializes the Fortran/C iteroperability types
 * @details This is called in the map_create_init_type routine following successful allocation 
 *          of memory. This should not be called directly by the user code. 
 * @param   init_type, Fortran/C interoperable type {@link MAP_InitInputType_t}
 * @see     map_create_init_type()
 */
void initialize_init_type_to_null(MAP_InitInputType_t* init_type);


/**
 * @brief   Initializes MAP internal initialization data
 * @details This is called in the map_create_init_type routine following successful allocation 
 *          of memory. This should not be called directly by the user code. The internal states
 *          are nullified and set to zero. 
 * @param   init_data, internal MAP initialization data {@link InitializationData}
 * @see     map_create_init_type()
 */
void initialize_init_data_to_null(InitializationData* init_data);


/**
 * @brief   Initialized all model data parameters to NULL or default values
 * @details Calls other functions to nullify and assign model data variables default values.          
 * @param   floater, the vessel strcture
 * @see     map_create_other_type()
 */
void initialize_domain_to_null(Domain* domain);


/**
 * @brief   Initializes MAP internal initialization data
 * @details This is called in the map_create_init_type routine following successful allocation 
 *          of memory. This should not be called directly by the user code. The internal states
 *          are nullified and set to zero. 
 * @param   options, model variables not related to solver items/options
 * @see     map_create_init_type()
 */
void initialize_model_option_defaults(DomainOptions* options);


/**
 * @brief   Sets inner-loop solver data 
 * @details Jacobian is 2x2 and residual vectors are 2x1. Sizes are fixed because only two
 *          algebraic equations are solved for each lines. 
 *          {@link free_outer_solve_data}
 * @param   inner, inner-loop solver data
 * @see     map_create_other_type()
 */
void initialize_inner_solve_data_defaults(InnerSolveAttributes* inner);


/**
 * @brief   Sets vessel parameters to NULL
 * @details This function nullifies Vessel-associated pointers. The Vessel structure is a
 *          convenient way to displace all fairleads connected to a vessel with a single 
 *          function call. Since the number of 'VESSEL' nodes are unknown at compile time,
 *          vessel arrays are sized at run-time. Deallocation occurs in {@link free_vessel}
 * @param   floater, the vessel strcture
 * @see     map_create_other_type()
 */
void initialize_vessel_to_null(Vessel* floater);


/**
 * @brief   Sets outer solver data (Jacobian, residual vector, ect) to NULL
 * @details Unlike the inner solver data strcture, Jacobian and residual vectors sizes are not 
 *          declared at compile time. These are resized depending on the number of 'connect'
 *          node declared in the problem. Outer solve data is deallocated in 
 *          {@link free_outer_solve_data}
 * @param   outer, outer-loop solver data
 * @see     map_create_other_type()
 */
void initialize_outer_solve_data_defaults(OuterSolveAttributes* outer);


/**
 * @brief   Sets model solver options corresponding to the MAP input file parameters 
 * @details Called in {@link map_init} to set the corresponding model (solver) options
 *          in MAP. This is different from line flags. Note that this function can 
 *          return MAP_WARNING and MAP_ERROR; these errors/warnings do not trip MAP_FATAL_33.
 *          Function must be modified to account for new options added at a future date. 
 *          Valid solver options input parameter include:
 *          <pre>
 *          help
 *          inner_ftol
 *          inner_gtol
 *          inner_xtol
 *          inner_max_its
 *          outer_max_its
 *          outer_tol
 *          outer_epsilon
 *          integration_dt (unsupported LM option)
 *          kb_default (unsupported LM option)
 *          cb_default (unsupported LM option)
 *          wave_kinematics (unsupported LM option)
 *          outer_bd
 *          outer_cd
 *          outer_fd
 *          pg_cooked (not compatible with krylov_accelerator option)
 *          krylov_accelerator (not compatible ith pg_cooked option)
 *          repeat
 *          ref_position
 *          </pre>
 * @param   domain, MAP interal data construct
 * @param   init_data, initialization input strings
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see     map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE set_model_options_list(Domain* data, InitializationData* init, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Initializes "HELP" flag model options
 * @details Called by {@link set_model_option_list} to set the boolean for printing help
 *          contents to stdout. Prefixing the line with a space (" ") treats the line 
 *          as a comment.
 *          MAP input file syntax:
 *          <pre>
 *          help
 *          </pre>
 * @param   list, a character array structure
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_help_flag(bstring list);


/**
 * @brief   Initializes "INNER_FTOL" tolerance for minpack
 * @details Called by {@link set_model_option_list} to set the function tolerance. 
 *          From MinPack documentation: "Termination occurs when both the actual 
 *          and predicted relative reductions in the sum of squares are at most FTOL. 
 *          Therefore, FTOL measures the relative error desired in the sum of squares."
 *          Prefixing the line with a space (" ") treats the line as a comment.
 *          MAP input file syntax:
 *          <pre>
 *          inner_ftol <float>
 *          </pre>
 * @param   list, a character array structure
 * @param   ftol, pointer to the FTOL paramter for minpack
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_inner_f_tol_flag(struct bstrList* list, double* ftol);


/**
 * @brief   Initializes "INNER_GTOL" tolerance for minpack
 * @details Called by {@link set_model_option_list} to set the function tolerance. 
 *          From MinPack documentation: "Termination occurs when the cosine of 
 *          the angle between FVEC and any column of the Jacobian is at most GTOL 
 *          in absolute value. Therefore, GTOL measures the orthogonality desired 
 *          between the function vector and the columns of the Jacobian." Prefixing
 *          the line with a space (" ") treats the line as a comment.
 *          MAP input file syntax:
 *          <pre>
 *          inner_gtol <float>
 *          </pre>
 * @param   list, a character array structure
 * @param   gtol, pointer to the GTOL paramter for minpack
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_inner_g_tol_flag(struct bstrList* list, double* gtol);


/**
 * @brief   Initializes "INNER_XTOL" tolerance for minpack
 * @details Called by {@link set_model_option_list} to set the function tolerance. 
 *          From MinPack documentation: "Termination occurs when the relative 
 *          error between two consecutive iterates is at most XTOL.  Therefore, 
 *          XTOL measures the relative error desired in the approximate solution.
 *          Prefixing the line with a space (" ") treats the line as a comment.
 *          MAP input file syntax:
 *          <pre>
 *          inner_xtol <float>
 *          </pre>
 * @param   list, a character array structure
 * @param   xtol, pointer to the XTOL paramter for minpack
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_inner_x_tol_flag(struct bstrList* list, double* xtol);


/**
 * @brief   Initializes "INNER_MAX_ITS" maximimum iteration count for inner-loop solver
 * @details Called by {@link set_model_option_list} to set the maximum number of 
 *          iterations for the MinPack routine. The solver terminate prematurely
 *          when this number is reached.
 *          MAP input file syntax:
 *          <pre>
 *          inner_max_its <float>
 *          </pre>
 * @param   list, a character array structure
 * @param   max_its, pointer to the maximum iteration count integer
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_inner_max_its_flag(struct bstrList* list, int* max_its);


/**
 * @brief   Initializes "OUTER_MAX_ITS" maximimum iteration count for outer-loop solver
 * @details Called by {@link set_model_option_list} to set the maximum number of 
 *          iterations for the MinPack routine. The solver terminate prematurely
 *          when this number is reaches
 *          MAP input file syntax:
 *          <pre>
 *          outer_max_its <float>
 *          </pre>
 * @param   list, a character array structure
 * @param   max_its, pointer to the maximum iteration count integer
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_outer_max_its_flag(struct bstrList* list, int* max_its);


/**
 * @brief   Initializes "OUTER_TOL" tolerance for custom outer-loopNewton-Raphson 
 *          solver
 * @details Called by {@link set_model_option_list} to set the function tolerance. 
 *          The convergence criteria required the norm of the function error to be 
 *          less than OUTER_TOL to Prefixing the line with a space (" ") treats the 
 *          line as a comment.
 *          MAP input file syntax:
 *          <pre>
 *          outer_tol <float>
 *          </pre>
 * @param   list, a character array structure
 * @param   outer_tol, pointer to the tolerance paramter
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_outer_tol_flag(struct bstrList* list, double* outer_tol);


/**
 * @brief   Initializes "OUTER_EPSILON" finite-difference increment
 * @details Called by {@link set_model_option_list} to set finite-difference tolerance
 *          for the outer-loop solver.
 *          MAP input file syntax:
 *          <pre>
 *          outer_epsilon <float>
 *          </pre>
 * @param   list, a character array structure
 * @param   epsilon, value by which the states are displaced
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_outer_epsilon_flag(struct bstrList* list, double* epsilon);


/**
 * @brief   Sets FD type to BACKWARD-DIFFERENCE 
 * @details Called by {@link set_model_option_list} to set the outer-loop finite-
 *          difference scheme to backward-difference. 
 *          $\frac{\mathbf{f}(x+\epsilon)-\mathbf{f}(x)}{\epsilon}$
 *          MAP input file syntax:
 *          <pre>
 *          outer_bd
 *          </pre>
 * @param   list, a character array structure
 * @param   bd, finite difference type
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_outer_bd_flag(struct bstrList* list, FdType* bd);

/**
 * @brief   Sets FdType type to CENTRAL-DIFFERENCE 
 * @details Called by {@link set_model_option_list} to set the outer-loop finite-
 *          difference scheme to central-difference. 
 *          $\frac{\mathbf{f}(x+\epsilon)-\mathbf{f}(x-\epsilon)}{2\epsilon}$
 *          MAP input file syntax:
 *          <pre>
 *          outer_bd
 *          </pre>
 * @param   list, a character array structure
 * @param   cd, finite difference type
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_outer_cd_flag(struct bstrList* list, FdType* cd);


/**
 * @brief   Sets FdType type to FORWARD-DIFFERENCE 
 * @details Called by {@link set_model_option_list} to set the outer-loop finite-
 *          difference scheme to forward-difference. 
 *          $\frac{\mathbf{f}(x+\epsilon)-\mathbf{f}(x)}{\epsilon}$
 *          MAP input file syntax:
 *          <pre>
 *          outer_fd
 *          </pre>
 * @param   list, a character array structure
 * @param   fd, finite difference type
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_outer_fd_flag(struct bstrList* list, FdType* fd);


/**
 * @brief   Sets DS and DSM preconditioning coefficients for outer-loop solver
 * @details Called by {@link set_model_option_list} to set preconditioner coefficients
 *          placed on the outer-loop Jacobian. With the cosntriant state update being:
 *          $\mathbf{D}_{i,i}=\frac{\mathbf{DS}}{n^{3/2}}+\mathbf{DSM}$
 *          then the preconditioner is:
 *          $\mathbf{z}_{n+1}=\mathbf{z}_{n}+\left(\mathbf{J}+\mathbf{D}\right)^{-1}\mathbf{r}$
 *          This is essentially a Jacobi preconditioner.
 *          MAP input file syntax:
 *          <pre>
 *          outer_fd
 *          </pre>
 * @article {peyrot1978,
 *           title={Analysis of flexible transmission lines},
 *           author={Peyrot, Alain H and Goulois, Alain M},
 *           journal={Journal of the Structural Division},
 *           volume={104},
 *           number={5},
 *           pages={763--779},
 *           year={1978},
 *           publisher={ASCE}
 *          }
 * @param   list, a character array structure
 * @param   solver, outer solvr attributes
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_pg_cooked_flag(struct bstrList* list, OuterSolveAttributes* solver);


/**
 * @brief   Sets flag to use Powell's method to solve the outer loop nonlinear equations
 *          MAP input file syntax:
 *          <pre>
 *          POWELL
 *          </pre>
 * @param   list, a character array structure
 * @param   solver, outer solvr attributes
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_powell_flag(struct bstrList* list, OuterSolveAttributes* solver);


/**
 * @brief   Sets integration time step for the LM model. Not supported yet!
 * @details Called by {@link set_model_option_list} to set the integrationt time step
 *          MAP input file syntax:
 *          <pre>
 *          integration_dt <float>
 *          </pre>
 * @param   list, a character array structure
 * @param   dt, integration time step
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_integration_dt_flag(struct bstrList* list, double* dt);


/**
 * @brief   Sets seabed stiffness factor for the LM model. Not supported yet!
 * @details Called by {@link set_model_option_list} to set the seabed stiffness
 *          MAP input file syntax:
 *          <pre>
 *          kb_default <float>
 *          </pre>
 * @param   list, a character array structure
 * @param   kb, seabed stiffness in [N/m]
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_kb_default_flag(struct bstrList* list, double* kb);


/**
 * @brief   Sets the boolean and maximum iterations for kyrlov accelerator
 * @details Called by {@link set_model_option_list} to choose the krylov
 *          accelerator option as the MSQS outer-loop option. Not compatible
 *          with pg_cooked option.
 *          MAP input file syntax:
 *          <pre>
 *          krylov_accelerator <float>
 *          </pre>
 * @param   list, a character array structure
 * @param   max_its, maximum numnber of krylov iterations 
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_krylov_accelerator_flag(struct bstrList* list, OuterSolveAttributes* solver);


/**
 * @brief   Sets seabed damping factor for the LM model. Not supported yet!
 * @details Called by {@link set_model_option_list} to set the seabed stiffness
 *          MAP input file syntax:
 *          <pre>
 *          cb_default <float>
 *          </pre>
 * @param   list, a character array structure
 * @param   kb, seabed damping in [N-s/m]
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_cb_default_flag(struct bstrList* list, double* cb);


/**
 * @brief   Compute the wave kinematics in MAP for the LM model. Not supported yet!
 * @details Called by {@link set_model_option_list} to enable calculation of the 
 *          local oribital fluid velocity and acceleration. This is fed as input
 *          to compute the added-mass and fluid-drag force. 
 *          MAP input file syntax:
 *          <pre>
 *          wave_kinematics
 *          </pre>
 * @param   list, a character array structure
 * @param   wave, true if wave kinematics are computed in MAP
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_wave_kinematics_flag(struct bstrList* list, bool* wave);


/**
 * @brief   Use the LM model. 
 * @details Called by {@link set_model_option_list} to set the mooring model to 
 *          the lumped-mass formulation. The MAP input file syntax:
 *          <pre>
 *          lm_model
 *          </pre>
 * @param   list, a character array structure
 * @param   lm, true if wave kinematics are computed in MAP
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_lm_model_flag(struct bstrList* list, bool* lm);


/**
 * @brief   Initializes "REPEAT" flag model options
 * @details Called by {@link set_model_option_list} to set the number of line and node 
 *          repeats. Prefixing the line with a space (" ") treats the line as a comment.
 *          MAP input file syntax:
 *          <pre>
 *          repeat <float_1> <float_2> ... <float_n>
 *          </pre>
 * @param   list, an array of bstring
 * @param   options, model options
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_repeat_flag(struct bstrList* list, DomainOptions* options);


/**
 * @brief   Set the model reference position 
 * @details Called by {@link set_model_option_list} to set the reference position. The
 *          reference position is where the linearized stiffness matrix is comuted. 
 *          <pre>
 *          re_position <float_1> <float_2> <float_3>
 *          </pre>
 * @param   list, an array of bstring
 * @param   ref_position, 3D reference position point
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_ref_position_flag(struct bstrList* list, Point* ref_position);


/**
 * @brief   Checks to make sure option is valid
 * @details Called by {@link set_model_option_list} to let users know if they have an 
 *          invalid option in te MAP input file. 
 * @param   list, an array of bstring
 * @see     set_model_options_list(), map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_uncaught_flag(struct bstrList* list);


/**
 * @brief   Sets cable properties
 * @details Called in {@link map_init} to set the cable properties. The specific properties
 *          do not necessarily need to be used. 
 *          <pre>
 *          label - cable name and identifier
 *          diam - cable diameter
 *          mass density (in air)
 *          EA - axial stiffness
 *          cb - seabed friction coefficient
 *          internal structural damping (unsupported LM option)
 *          cross-flow added mass coefficient (unsupported LM option)
 *          cross-float damping coefficient (unsupported LM option)
 *          tangen drag coefficient (unsupported LM option)
 *          </pre>
 * @param   domain, MAP interal data construct
 * @param   init_data, initialization input strings
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see     map_init()
 * @return  MAP error code
 */
MAP_ERROR_CODE set_cable_library_list( Domain* domain, InitializationData* init_data, char* map_msg, MAP_ERROR_CODE* ierr );


/**
 * @brief   Set line diameter at initialization 
 * @param   word, word to conver to float
 * @param   library_ptr, pointer to a cable library
 * @see     set_cable_library_list()
 * @return  MAP error code
 */
MAP_ERROR_CODE set_library_diameter(bstring word, CableLibrary* library_ptr);


/**
 * @brief   Set line mass density in air at initialization 
 * @param   word, word to conver to float
 * @param   library_ptr, pointer to a cable library
 * @see     set_cable_library_list()
 * @return  MAP error code
 */
MAP_ERROR_CODE set_library_mass_density(bstring word, CableLibrary* library_ptr);


/**
 * @brief   Set line axial stiffness at initialization 
 * @param   word, word to conver to float
 * @param   library_ptr, pointer to a cable library
 * @see     set_cable_library_list()
 * @return  MAP error code
 */
MAP_ERROR_CODE set_library_ea(bstring word, CableLibrary* library_ptr);


/**
 * @brief   Set line/seabed friction coefficient at initialization 
 * @param   word, word to conver to float
 * @param   library_ptr, pointer to a cable library
 * @see     set_cable_library_list()
 * @return  MAP error code
 */
MAP_ERROR_CODE set_library_cb(bstring word, CableLibrary* library_ptr);


/**
 * @brief   Set line structural damping (internal) coefficient at initialization 
 * @param   word, word to conver to float
 * @param   library_ptr, pointer to a cable library
 * @see     set_cable_library_list()
 * @return  MAP error code
 */
MAP_ERROR_CODE set_library_internal_damping(bstring word, CableLibrary* library_ptr);


/**
 * @brief   Set line cross-flow added mass coefficient initialization 
 * @param   word, word to conver to float
 * @param   library_ptr, pointer to a cable library
 * @see     set_cable_library_list()
 * @return  MAP error code
 */
MAP_ERROR_CODE set_library_added_mass_coefficient(bstring word, CableLibrary* library_ptr);


/**
 * @brief   Set line cross-flow drag coefficient coefficient initialization 
 * @param   word, word to conver to float
 * @param   library_ptr, pointer to a cable library
 * @see     set_cable_library_list()
 * @return  MAP error code
 */
MAP_ERROR_CODE set_library_cross_flow_drag_coefficient(bstring word, CableLibrary* library_ptr);


/**
 * @brief   Set line skin drag coefficient initialization 
 * @param   word, word to conver to float
 * @param   library_ptr, pointer to a cable library
 * @see     set_cable_library_list()
 * @return  MAP error code
 */
MAP_ERROR_CODE set_library_tangent_drag_coefficient(bstring word, CableLibrary* library_ptr);


/**
 * @brief   Reset line properties to null or zero
 * @details The only reason this is necessary is to destroy the label parameter
 *          in the cable library structure. If this parameter is not destroyed, 
 *          then an error occurs when the next cable library entry is read. A 
 *          segfault may occur when the bstrcpy() is called a second time on a
 *          non-nullified bstring. 
 * @param   word, word to conver to float
 * @param   library_ptr, pointer to a cable library
 * @see     set_cable_library_list()
 * @return  MAP error code
 */
MAP_ERROR_CODE reset_cable_library(CableLibrary* new_cable_library);


MAP_ERROR_CODE repeat_nodes(Domain* domain, InitializationData* init_data, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE set_node_list(const MAP_ParameterType_t* p_type,  MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, Domain* domain, struct bstrList* node_input_string, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE allocate_types_for_nodes(MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, Domain* domain, struct bstrList* node_input_string, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE reset_node(Node* node_ptr);
MAP_ERROR_CODE expand_node_number(const int n_line, bstring line); 
MAP_ERROR_CODE expand_node_type(const char* word, bstring line); 
MAP_ERROR_CODE expand_node_position_x(double* x, const char* word);
MAP_ERROR_CODE expand_node_position_y(double* y, const char* word);
MAP_ERROR_CODE expand_node_position_z(Vector* position, const double angle, const double x, const double y, const char* word, bstring line);
MAP_ERROR_CODE expand_node_mass(const char* word, bstring line); 
MAP_ERROR_CODE expand_node_buoyancy(const char* word, bstring line); 
MAP_ERROR_CODE expand_node_force_x(double* fx, const char* word); 
MAP_ERROR_CODE expand_node_force_y(double* fy, const char* word); 
MAP_ERROR_CODE expand_node_force_z(Vector* force, const double angle, const double fx, const double fy, const char* word, bstring line);


MAP_ERROR_CODE repeat_lines(Domain* domain, InitializationData* init, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE set_line_list(MAP_ConstraintStateType_t* z_type, Domain* domain, struct bstrList* line_input_string, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE reset_line(Line* line_ptr);
MAP_ERROR_CODE expand_line_number(const int n_line, bstring line); 
MAP_ERROR_CODE expand_line_property_name(const char* word, bstring line); 
MAP_ERROR_CODE expand_line_length(const char* word, bstring line); 
MAP_ERROR_CODE expand_line_anchor_number(const char* word, const int index, const int n, bstring line); 
MAP_ERROR_CODE expand_line_fairlead_number(const char* word, const int index, const int n, bstring line); 
MAP_ERROR_CODE expand_line_flag(const char* word, bstring line); 


MAP_ERROR_CODE set_vartype_ptr(const char* unit, bstring alias, const int num, VarTypePtr* type, bstring property);
MAP_ERROR_CODE set_vartype(const char* unit, bstring alias, const int num, VarType* type, bstring property );


/**
 * @brief   Compares two integers.
 * @details Used in comparing the size of a MAP state against iteration count. 
 *          Used in {@link set_node_list()}
 * @param   a, an integer 
 * @param   b, an integer
 * @return  MAP_FATAL if a != b; otherwise safe
 */
MAP_ERROR_CODE compare_integer_length(const int a, const int b);


// MAP_ERROR_CODE set_line_vartype(Line* line_ptr, const int i); @rm
MAP_ERROR_CODE associate_line_with_cable_property(Line* line_ptr, Domain* domain, const char* word, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE associate_line_with_anchor_node(Line* line_ptr, Domain* domain, const int line_num, const char* word, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE associate_line_with_fairlead_node(Line* line_ptr, Domain* domain, const int line_num, const char* word, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE set_line_option_flags(struct bstrList* word, int* index, Line* line_ptr, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   set the 'unit' and 'name' attributes for the node vartype's
 * @details Used in the initialization routine, part of {@link set_node_list()}
 * @param   node_ptr pointer to a newly allocated node data struct
 */
MAP_ERROR_CODE set_node_vartype(Node* node_ptr);


/**
 * @brief
 * @param domain internal state data {@link Domain}
 * @param io_type initialization output type, native C struct {@link MAP_InitOutputType_t}
 * @param map_msg error message
 * @param ierr error code
 */
MAP_ERROR_CODE set_output_list(Domain* domain, MAP_InitOutputType_t* io_type, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief This sets the pointers to NULL for the vessel object and gives it default properties. 
 * @todo: need to associate the node with inputs
 * @acceses: set_vartype_float( )
 * @calledby: mapcall_msqs_init( )
 */
MAP_ERROR_CODE set_vessel(Vessel* floater, const MAP_InputType_t* u_type, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * This sets the pointers to NULL for the vessel object and gives it default properties. Only 
 * to be used in the python glue code. 
 *
 * @todo: need to associate the node with inputs
 * @acceses: set_vartype_float( )
 * @calledby: mapcall_msqs_init( )
 */
MAP_ERROR_CODE set_vartype_float(const char* unit, const char* alias, const int num, VarType* type, const double value);


/**
 * @breif   Allocate outer solve data
 * @param   ns, outer solve attributes; preserves l, jac, and u
 * @param   size, matrix size
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see mapcall_msqs_init()
 * @see mapcall_msqs_end()
 * @see free_outer_solve_data() for data deallocation
 * @return  MAP error code
 */
MAP_ERROR_CODE allocate_outer_solve_data(OuterSolveAttributes* ns, const int size, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @breif   Initialized omega (weight per unit length) and cross-section area of a cable. 
 * @details Called by {@link map_init} to assign:
 *          $A=\pi*\frac{radius^{2}}{4}$
 *          and
 *          $\omega=g*(\mu-A*\rho_{seawater})$ 
 * @param   domain, mooring domain
 * @param   p_type, parameter type
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see mapcall_msqs_init()
 * @see mapcall_msqs_end()
 * @see free_outer_solve_data() for data deallocation
 * @return  MAP error code
 */
MAP_ERROR_CODE initialize_cable_library_variables(Domain* domain, MAP_ParameterType_t* p_type, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * Sets init data to NULL or -9999
 *
 * @acceses: none
 * @calledby: mapcall_msqs_init( )
 *
 *  OBTAINED FROM THE CMINPACK SOURCE. The following defintions of the inputs for lmder(...) are provided in the cminpack docmentation.  
 *  
 *  info = __cminpack_func__(lmder)( inner_function_evals, lineIter, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, 
 *                                   maxfev, diag, mode, factor, nprint, &nfev, &njev, ipvt, qtf, wa1, wa2, wa3, wa4);
 *     
 *      - ftol          : a nonnegative input variable. Termination occurs when both the actual and predicted relative reductions in the sum of squares are at most ftol. Therefore, ftol measures the relative error desired in the sum of squares.
 *      - xtol          : a nonnegative input variable. Termination occurs when the relative error between two consecutive iterates is at most xtol. Therefore, xtol measures the relative error desired in the approximate solution.
 *      - gtol          : a nonnegative input variable. Termination occurs when the cosine of the angle between fvec and any column of the Jacobian is at most gtol in absolute value. Therefore, gtol measures the orthogonality desired between the function vector and the columns of the Jacobian.
 *      - maxfev        : a positive integer input variable. Termination occurs when the number of calls to fcn is at least maxfev by the end of an iteration.
 *      - diag          : an array of length n. If mode = 1 (see below), diag is internally set. If mode = 2, diag must contain positive entries that serve as multiplicative scale factors for the variables.
 *      - mode          : an integer input variable. If mode = 1, the variables will be scaled internally. If mode = 2, the scaling is specified by the input diag. Other values of mode are equivalent to mode = 1.
 *      - factor        : a positive input variable used in determining the initial step bound. This bound is set to the product of factor and the euclidean norm of diag*x if the latter is nonzero, or else to factor itself. In most cases factor should lie in the interval (.1,100.). 100. is a generally recommended value.
 *      - nprint        : an integer input variable that enables controlled printing of iterates if it is positive. In this case, fcn is called with iflag = 0 at the beginning of the first iteration and every nprint iterations thereafter and immediately prior to return, with x and fvec available for printing. If nprint is not positive, no special calls of fcn with iflag = 0 are made.
 *      - info          : an integer output variable. If the user has terminated execution, info is set to the (negative) value of iflag. See description of fcn. Otherwise, info is set as follows.
 *      - info=0        : improper input parameters.
 *      - info=1        : both actual and predicted relative reductions in the sum of squares are at most ftol.
 *      - info=2        : relative error between two consecutive iterates is at most xtol.
 *      - info=3        : conditions for info = 1 and info = 2 both hold.
 *      - info=4        : the cosine of the angle between fvec and any column of the Jacobian is at most gtol in absolute value.
 *      - info=5        : number of calls to fcn has reached or exceeded maxfev.
 *      - info=6        : ftol is too small. No further reduction in the sum of squares is possible.
 *      - info=7        : xtol is too small. No further improvement in the approximate solution x is possible.
 *      - info=8        : gtol is too small. fvec is orthogonal to the columns of the Jacobian to machine precision.
 *      - nfev          : an integer output variable set to the number of calls to fcn with iflag = 1.
 *      - njev          : an integer output variable set to the number of calls to fcn with iflag = 2.
 *      - ipvt          : an integer output array of length n. ipvt defines a permutation matrix p such that jac*p = q*r, where jac is the final calculated Jacobian, q is orthogonal (not stored), and r is upper triangular with diagonal lines of nonincreasing magnitude. Column j of p is column ipvt(j) of the identity matrix.
 *      - qtf           : an output array of length n which contains the first n lines of the vector (q transpose)*fvec.
 *      - wa1, wa2, wa3 : are work arrays of length n.
 *      - wa4           : a work array of length m.  
 */
MAP_ERROR_CODE first_solve(Domain* domain, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * Check if a string can be converted to a valid numeric value
 *
 * 0 (MAP_SAFE) is reserved for 'safe' exit out of functions. 
 *
 * @param s     must point to a user-provided memory location
 * @return      0 for success, 3 for failure
 *
 * @see MAP_ERROR_CODE
 */
MAP_ERROR_CODE is_numeric(const char* string);


/**
 * @breif   Write the map summary file
 * @details Is called once in {@link map_init} to create the summary file post-solve.
 *          Even if the solution is not solved, the summary file will be written. 
 * @param   init_type, initialization type
 * @param   p_type, parameter type
 * @param   y_type, output type
 * @param   other_type, oter state type
 * @param   domain, mooring domain
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see mapcall_msqs_init()
 */
void log_initialization_information(MAP_InitInputType_t* init_type, MAP_ParameterType_t* p_type, MAP_OutputType_t* y_type, MAP_OtherStateType_t* other_type, Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr);


MAP_ERROR_CODE associate_vartype_ptr(VarTypePtr* type, double* arr, int index);
void copy_target_string(char* target, unsigned char* source);
MAP_ERROR_CODE initialize_external_applied_force(char* unit, char* alias, const int num, VarType* type, char const* property);
MAP_ERROR_CODE initialize_node_sum_force_ptr(char* unit, char* alias, const int num, VarTypePtr* type);
MAP_ERROR_CODE map_get_version(MAP_InitOutputType_t* io_type);
MAP_ERROR_CODE print_help_to_screen();
size_t vartype_meter(const void* el);
size_t vartype_ptr_meter(const void *el);
const char* remove_first_character(const char* string);
void print_machine_name_to_screen();

/**
 */
MAP_ERROR_CODE associate_constraint_states(Domain* domain, MAP_ConstraintStateType_t* z_type);


#endif /* _INITIALIZATION_H */
