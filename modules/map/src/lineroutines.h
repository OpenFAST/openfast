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


#ifndef _LINEROUTINES_H
#define _LINEROUTINES_H


#include "map.h"


/**
 * success = incremenet_x_dof_by_delta(uType, -epsilon);
 */
MAP_ERROR_CODE increment_dof_by_delta(double* u_type, const double delta, const int size);


/**
 * success = restore_original_displacement(uType->x, xOriginal, N);
 */
MAP_ERROR_CODE restore_original_displacement(double* u_ype, const double* initial_value, const int size);


/**
 * success = reset_force_to_zero(yType, N);
 */
MAP_ERROR_CODE reset_force_to_zero(double* fx, double* fy, double* fz, double* mx, double* my, double* mz, const int size);


/**
 * success = set_force_plus(yType->Fx, fx, N);
 */
MAP_ERROR_CODE set_force_minus(const double* in_force, double* force, const int size);


/**
 * success = set_force_plus(yType->Fx, fx, N);
 */
MAP_ERROR_CODE set_force_plus(const double* in_force, double* force, const int size);


/**
 *
 */
MAP_ERROR_CODE update_outer_loop_inputs(double* input, MAP_ConstraintStateType_t* z_type,  const int size, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 *
 */
MAP_ERROR_CODE update_outer_loop_residuals(double* residual, MAP_OtherStateType_t* other_type, const int size, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * success = set_moment_minus(yType, vessel, mx, my, mz, N)
 */
MAP_ERROR_CODE set_moment_minus(const MAP_OutputType_t* y_type, const Vessel* vessel, double* mx, double* my, double* mz, const int size);


MAP_ERROR_CODE set_moment_plus(const MAP_OutputType_t* y_type, const Vessel* vessel, double* mx, double* my, double* mz, const int size);


MAP_ERROR_CODE set_moment_minus_phi(const MAP_InputType_t* u_type, const MAP_OutputType_t* y_type, const Vessel* vessel, double* mx, double* my, double* mz, const double epsilon, const int size);
MAP_ERROR_CODE set_moment_plus_phi(const MAP_InputType_t* u_type, const MAP_OutputType_t* y_type, const Vessel* vessel, double* mx, double* my, double* mz, const double epsilon, const int size);
MAP_ERROR_CODE set_moment_minus_the(const MAP_InputType_t* u_type, const MAP_OutputType_t* y_type, const Vessel* vessel, double* mx, double* my, double* mz, const double epsilon, const int size);
MAP_ERROR_CODE set_moment_plus_the(const MAP_InputType_t* u_type, const MAP_OutputType_t* y_type, const Vessel* vessel, double* mx, double* my, double* mz, const double epsilon, const int size);
MAP_ERROR_CODE set_moment_minus_psi(const MAP_InputType_t* u_type, const MAP_OutputType_t* y_type, const Vessel* vessel, double* mx, double* my, double* mz, const double epsilon, const int size);
MAP_ERROR_CODE set_moment_plus_psi(const MAP_InputType_t* u_type, const MAP_OutputType_t* y_type, const Vessel* vessel, double* mx, double* my, double* mz, const double epsilon, const int size);


MAP_ERROR_CODE increment_phi_dof_by_delta(MAP_InputType_t* u_type, const Vessel* vessel, const double delta, const int size);


MAP_ERROR_CODE increment_the_dof_by_delta(MAP_InputType_t* u_type, const Vessel* vessel, const double delta, const int size);


MAP_ERROR_CODE increment_psi_dof_by_delta(MAP_InputType_t* u_type, const Vessel* vessel, const double delta, const int size);


MAP_ERROR_CODE f_op_sequence(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, int size, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_x_sequence(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_pos, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_y_sequence(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_pos, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_z_sequence(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_pos, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_phi_sequence(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_x, const double* original_y, const double* original_z, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_the_sequence(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_x, const double* original_y, const double* original_z, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_psi_sequence(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, Fd* force, const double epsilon, const int size, const double* original_x, const double* original_y, const double* original_z, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE calculate_stiffness(double* K, Fd* force, const double delta, const int size);
MAP_ERROR_CODE calculate_sumforce (double* F, Fd* force, const int size);

/**
 * sets cable excursions (l and h) and reference frame psi rotation
 */
MAP_ERROR_CODE set_line_variables_pre_solve(Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * set: alpha, alpha at anchor, H, V and T at anchor, T at fairlead
 */
MAP_ERROR_CODE set_line_variables_post_solve(Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr);

double set_vertical_excursion(Line* line);
double set_horizontal_excursion(Line* line);


/**
 *  'psi' is the angle between the line x-axis (local frame) and X-axis (global frame). This essentially produces this rotation matrix:
 *  
 *  \mathbf{R}(\psi) = \begin{bmatrix}
 *                     \cos\psi & -\sin\psi & 0 \\ 
 *                     \sin\psi &  \cos\psi & 0 \\ 
 *                            0 &         0 & 1
 *                     \end{bmatrix}
 *  
 *       1) first find psi - the angle of rotation between the line frame and the global reference frame
 *       2) r_j = fairlead displacement
 *       3) r_i = anchor displacement
 *       4) cos^{-1}( dot( (r_j-r_i) , (u_i) ) / ( norm(r_j-r_i) ) )
 */
MAP_ERROR_CODE set_psi(Line* line, char* map_msg, MAP_ERROR_CODE* ierr);


MAP_ERROR_CODE reset_node_force_to_zero(Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * Initialized omega (weight per unit length) and cross-section area of a cable. The formula is 
 *
 *   A=\pi*\frac{radius^{2}}{4}
 *   \omega=g*(\mu-A*\rho_{seawater})
 *
 * @acceses: none
 * @calledby: mapcall_msqs_init( )
 */                                                   
MAP_ERROR_CODE set_line_initial_guess(Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr);


MAP_ERROR_CODE solve_linear_spring_cable(Line* line, char* map_msg, MAP_ERROR_CODE* ierr);

/**
 * MAP_InputType_t* uType,
 * MAP_ConstraintStateType_t* zType,
 * MAP_OtherStateType_t* otherType,
 * MAP_OutputType_t* yType,
 */
MAP_ERROR_CODE node_solve_sequence(Domain* domain, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, const float time, char* map_msg, MAP_ERROR_CODE* ierr);

MAP_ERROR_CODE newton_solve_sequence(Domain* domain, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, const float time, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE krylov_solve_sequence(Domain* domain, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, const float time, char* map_msg, MAP_ERROR_CODE* ierr);

MAP_ERROR_CODE line_solve_sequence(Domain* domain, MAP_ParameterType_t* p_type, const float time, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE solve_line(Domain* domain, const float time, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Determine maximum line length before double-backing occurs
 * @details Called by {@link solve_line} to make sure the line is not double-
 *          backing on itself. In the extreme case, the cable resembles an 
 *          'L' shape. When the rule is violated, the upper tipe of the 'L'
 *          bends to the left. A violation of the algebrais equations happens. 
 *          You can disregard double-backing error checking by raising the 
 *          'omit_contact' element flag.
 *          Double backing is present when:
 *          $L_{u} >= l - \frac{EA}{\omega} + \sqrt{(\frac{EA}{W}^2) + 2*h*\frac{EA}{W}};
 *          This is an adaptation from the following Fortran code:
 *          <pre>
 *          ELSEIF ( W  >  0.0_DbKi )  THEN   ! .TRUE. when the line will sink in fluid
 *             LMax      = XF - EA/W + SQRT( (EA/W)*(EA/W) + 2.0_DbKi*ZF*EA/W )  ! Compute the maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead)
 *             IF ( ( L  >=  LMax   ) .AND. ( CB >= 0.0_DbKi ) )  &  ! .TRUE. if the line is as long or longer than its maximum possible value with seabed interaction
 *                CALL ProgAbort ( ' Unstretched mooring line length too large. '// &
 *                                 ' Routine Catenary() cannot solve quasi-static mooring line solution.' )
 *          ENDIF
 *          </pre>
 * @param   line, the line structure to be checked
 * @param   contact_flag, flag to ignore the double-backing check
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see     solve_line()
 * @return  MAP error code
 */
MAP_ERROR_CODE check_maximum_line_length(Line* line, const bool contact_flag, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * call immediately after set_line_variables_post_solve(); this added H and V
 */
MAP_ERROR_CODE calculate_node_sum_force(Domain* domain, MAP_ParameterType_t* p_type);



/**
 * increment sum force value by (f) if node is fairlead; (-f) is node is anchor. 
 */
void add_to_sum_fx(Node* node, const double fx);


/**
 * increment sum force value by (f) if node is fairlead; (-f) is node is anchor. 
 */
void add_to_sum_fy(Node* node, const double fy);


/**
 * increment sum force value by (f) if node is fairlead; (-f) is node is anchor. 
 */
void add_to_sum_fz(Node* node, const double fz);


#endif /* _LINEROUTINES_H */
