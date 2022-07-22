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


#ifndef _MAPAPI_H
#define _MAPAPI_H


/**
 * @brief   Initalizes all MAP base types (including some internal state)
 * @details The idea is to set variables to zero and null to prevent seg-faults in the case of 
 *          early program termination before initialization (MAP_Init) is fully complete. 
 *          {@link MAP_Init}
 * @param   el, opaque object used in simclist
 * @see     map_init()
 * @return  Size of CableLibrary structure
 */
MAP_EXTERNCALL void map_initialize_msqs_base(MAP_InputType_t* u_type, MAP_ParameterType_t* p_type, MAP_ContinuousStateType_t* x_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, MAP_InitOutputType_t* io_type);


MAP_EXTERNCALL void set_init_to_null(MAP_InitInputType_t* init_type, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL void map_add_cable_library_input_text(MAP_InitInputType_t* init_type);
MAP_EXTERNCALL void map_add_node_input_text(MAP_InitInputType_t* init_type);
MAP_EXTERNCALL void map_add_line_input_text(MAP_InitInputType_t* init_type);
MAP_EXTERNCALL void map_add_options_input_text(MAP_InitInputType_t* init_type);
MAP_EXTERNCALL double* map_plot_x_array(MAP_OtherStateType_t* other_type, int i, int num_points, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL double* map_plot_y_array(MAP_OtherStateType_t* other_type, int i, int num_points, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL double* map_plot_z_array(MAP_OtherStateType_t* other_type, int i, int num_points, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL void map_plot_array_free(double* array) ;
MAP_EXTERNCALL double map_residual_function_length(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL double map_residual_function_height(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL double map_jacobian_dxdh(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL double map_jacobian_dxdv(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL double map_jacobian_dzdh(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL double map_jacobian_dzdv(MAP_OtherStateType_t* other_type, int i, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_EXTERNCALL int map_size_lines(MAP_OtherStateType_t* other_type, MAP_ERROR_CODE* ierr, char* map_msg);


/**
 * @brief     Deallocates the memory space for the init structure. Should be called immediately after map_init()
 * @param     init MAP-native initialization data structure. This is distinct from the FAST-framework data structure 
 * @param     map_msg error string
 * @param     ierr error code
 * @return    MAP_SAFE if it completes successfully
 * @see       {@link map_init()}
 */
MAP_EXTERNCALL int free_init_data (InitializationData* init, char* map_msg, MAP_ERROR_CODE* ierr);

/**
 * @brief Set the water depth. Should be called before {@link map_init()}
 * @param p_type paramter type, native C struct {@link MAP_ParameterType_t}
 * @param depth water depth [m]
 *
 * Example Fortran usage:
 * @code
 * INTERFACE                                                                
 *    SUBROUTINE mapextern_set_depth(interf, fc_val) BIND(C,name='map_set_sea_depth')  
 *      IMPORT                            
 *      IMPLICIT NONE                     
 *      TYPE(MAP_ParameterType_C) interf
 *      REAL(C_DOUBLE), VALUE :: fc_val  
 *    END SUBROUTINE mapextern_set_depth    
 * END INTERFACE                      
 *
 *   ...
 *
 * ! access the function using this subroutine call: 
 * CALL mapextern_set_depth(p%C_obj, depth)
 * @endcode
 */
MAP_EXTERNCALL void map_set_sea_depth(MAP_ParameterType_t* p_type, const double depth);


/**
 * @brief Set the water density. Should be called before {@link map_init()}
 * @param p_type paramter type, native C struct {@link MAP_ParameterType_t}
 * @param rho water density [kg/m^3]
 *
 * Example Fortran usage:
 * @code
 * INTERFACE                                                                         
 *    SUBROUTINE mapextern_set_density(interf, fc_val) BIND(C,name='map_set_sea_density') 
 *      IMPORT                                            
 *      IMPLICIT NONE                                     
 *      TYPE(MAP_ParameterType_C) interf                
 *      REAL(C_DOUBLE), VALUE :: fc_val                      
 *    END SUBROUTINE mapextern_set_density                      
 * END INTERFACE                                          
 *
 *   ...
 *
 * ! access the function using this subroutine call: 
 * CALL mapextern_set_density(p%C_obj, rho)
 * @endcode
 */
MAP_EXTERNCALL void map_set_sea_density(MAP_ParameterType_t* p_type, const double rho);


/**
 * @brief Set the gravitational constant. Should be called before {@link map_init()}
 * @param p_type paramter type, native C struct {@link MAP_ParameterType_t}
 * @param grtavity gravitational acceleration [m/s^2]
 *
 * Example Fortran usage:
 * @code
 * INTERFACE                                                               
 *    SUBROUTINE MAP_map_set_gravity(interf, fc_val) BIND(C,name='map_set_gravity')
 *      IMPORT                                                         
 *      IMPLICIT NONE                                                  
 *      TYPE(MAP_ParameterType_C) interf                               
 *      REAL(C_DOUBLE), VALUE :: fc_val                                
 *    END SUBROUTINE MAP_map_set_gravity                                   
 * END INTERFACE                                                       
 *
 *   ...
 *
 * ! access the function using this subroutine call: 
 * CALL mapextern_map_set_gravity(p%C_obj, g)
 * @endcode
 */
MAP_EXTERNCALL void map_set_gravity(MAP_ParameterType_t* p_type, const double gravity);


/**
 * @brief   Returns vertical and horizontal fairlead force along line plane
 * @details 
 * @param   H, reference to horizontal fairlead force magnitude
 * @param   V, reference to vertical fairlead force magnitude
 * @param   other_type, other state type fortran derived
 * @param   index, line number we are requesting the data for
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see     
 */
MAP_EXTERNCALL void map_get_fairlead_force_2d(double* H, double* V, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Returns X, Y, and Z fairlead force in global reference frame
 * @details
 * @param   fx, reference to horizontal X fairlead force in global frame
 * @param   fy, reference to horizontal Y fairlead force in global frame
 * @param   fz, reference to vertical Z fairlead force in global frame
 * @param   other_type, other state type fortran derived
 * @param   index, line number we are requesting the data for
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see     
 */
MAP_EXTERNCALL void map_get_fairlead_force_3d(double* fx, double* fy, double* fz, MAP_OtherStateType_t* other_type, int index, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * call this in python: offset_vessel().argtypes = [MapData_Type, MapInput_Type, c_char_p, POINTER(c_int)]        
 * angles are radians
 * 
 *     lib.map_offset_vessel.argtypes = [MapData_Type, MapInput_Type, c_double, c_double, c_double, c_double, c_double, c_double, c_char_p, POINTER(c_int)]        
 *
 * Angles are in degrees. This routine converts angles from deg to rad
 * Per Fossen (Fossen, Thor I. Guidance and control of ocean vehicles. Vol. 199. No. 4. New York: Wiley, 1994), this 
 * matrix converts vector from the body frame to the global reference frame:
 *
 * R = | cos(ψ)*cos(θ)    cos(ψ)*sin(θ)*sin(φ) − sin(ψ)*cos(φ)       cos(ψ)*sin(θ)*cos(φ) + sin(ψ)*sin(φ)  |
 *     | sin(ψ)*cos(θ)    sin(ψ)*sin(θ)*sin(φ) + cos(ψ)*cos(φ)       sin(ψ)*sin(θ)*cos(φ) − cos(ψ)*sin(φ)  |
 *     |   −sin(θ)                   cos(θ)*sin(φ)                                cos(θ)*cos(φ)            |
 *
 * We need to invoke this procedure to move the vessel nodes with body rotations factored:
 *
 * u_type = x + [R]*r       
 *   ▲     ▲       ▲
 *   |     |       |
 * global  |       |
 *      global     |
 *               local
 */
MAP_EXTERNCALL void map_offset_vessel(MAP_OtherStateType_t* other_type, MAP_InputType_t* u_type, double x, double y, double z, double phi, double the, double psi, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * lib.linearize_matrix.argtypes = [MapInput_Type, MapData_Type, MapOutnput_Type, c_double, c_char_p, POINTER(c_int)]        
 */
MAP_EXTERNCALL double** map_linearize_matrix(MAP_InputType_t* u_type, MAP_ParameterType_t* p_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, double epsilon, MAP_ERROR_CODE* ierr, char* map_msg);


/**
 * lib.py_free_linearize_matrix.argtypes = [POINTER(POINTER(c_double))]
 */
MAP_EXTERNCALL void map_free_linearize_matrix(double** array);

/**
 * lib.map_f_op.argtypes = [MapInput_Type, MapData_Type, MapOutnput_Type, c_double, c_char_p, POINTER(c_int)]        
 */
MAP_EXTERNCALL double* map_f_op(MAP_InputType_t* u_type, MAP_ParameterType_t* p_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, MAP_ConstraintStateType_t* z_type, MAP_ERROR_CODE* ierr, char* map_msg);

/**
 * lib.map_free_f_op.argtypes = [POINTER(c_double)]
 */
MAP_EXTERNCALL void map_free_f_op(double* array);


/**
 * @brief     Initializes the MAP model and allocates memory. Inconsistencies with the input file is reported here. 
 *            This should only be called once. 
 * @param     init_type initialization type, F2C FAST-native derived type
 * @param     u_type input type, F2C FAST-native derived type. 
 * @param     p_type parameter type, F2C FAST-native derived type
 * @param     x_type continuous-state type, F2C FAST-native derived type
 * @param     none void type. It exists because MAP does not use discrete types (it would fill this place otherwise)
 * @param     z_type constraint-state type, F2C FAST-native derived tpe
 * @param     other_type other-state type, F2C FAST-native derived type. This is core structure that houses all internal MAP states that cannot be mapped to Fortran
 * @param     y_type output type, F2C FAST-native derived type
 * @param     ioType init-output types, F2C FAST-native derived type
 * @param     ierr error code
 * @param     map_msg error string
 */
MAP_EXTERNCALL void map_init(MAP_InitInputType_t* init_type, 
                             MAP_InputType_t* u_type,
                             MAP_ParameterType_t* p_type,
                             MAP_ContinuousStateType_t* x_type,
                             MAP_DiscreteStateType_t* xd_type,
                             MAP_ConstraintStateType_t* z_type,
                             MAP_OtherStateType_t* other_type,
                             MAP_OutputType_t* y_type,
                             MAP_InitOutputType_t* ioType,
                             MAP_ERROR_CODE* ierr,
                             char* map_msg);


/**
 * @brief     Solves the statics problem for the MSQS system and should be called at each time step or vessel displacement. 
 * @details   Can be called multiple times, but must be called between {@link map_init()} and {@link map_end()} 
 *            If the reference to u_type changes, then we have to update the location MAP internal states are pointing 
 *            to. This is accomplished in the following code. The issue here is when this is called in Fortran:
 * 
 *                CALL MAP_CopyInput(u(1), u_interp, MESH_NEWCOPY, ErrStat, ErrMsg)      
 *
 *            u_interp is passed into into the argument for map_update_states(); however, the internal states are not
 *            pointing to data in u_interp. We address this below. Note that the initial reference for point_iter is set
 *            in {@link set_node_list()}    
 * @param     t current (global) time
 * @param     interval coupling interval
 * @param     u_type input type, F2C FAST-native derived type. 
 * @param     p_type parameter type, F2C FAST-native derived type
 * @param     x_type continuous-state type, F2C FAST-native derived type
 * @param     none void type. It exists because MAP does not use discrete types (it would fill this place otherwise)
 * @param     z_type constraint-state type, F2C FAST-native derived tpe
 * @param     other_type other-state type, F2C FAST-native derived type. This is core structure that houses all internal MAP states that cannot be mapped to Fortran
 * @param     ierr error code
 * @param     map_msg error string
 * @see       {@link map_calc_output()}
 */
MAP_EXTERNCALL void map_update_states(float t,
                                      int interval,
                                      MAP_InputType_t* u_type,
                                      MAP_ParameterType_t* p_type,
                                      MAP_ContinuousStateType_t* x_type,
                                      MAP_DiscreteStateType_t* xd_type,
                                      MAP_ConstraintStateType_t* z_type,
                                      MAP_OtherStateType_t* other_type,
                                      MAP_ERROR_CODE* ierr,
                                      char* map_msg);


/**
 * @brief     Retrieves the values after the statics problem is solved. This function should be called immediately after map_update_states. 
 *            Can be called multiple times, but must be called between {@link map_init()} and {@link map_end()}
 * @param     t current (global) time
 * @param     u_type input type, F2C FAST-native derived type. 
 * @param     p_type parameter type, F2C FAST-native derived type
 * @param     x_type continuous-state type, F2C FAST-native derived type
 * @param     none void type. It exists because MAP does not use discrete types (it would fill this place otherwise)
 * @param     z_type constraint-state type, F2C FAST-native derived tpe
 * @param     other_type other-state type, F2C FAST-native derived type. This is core structure that houses all internal MAP states that cannot be mapped to Fortran
 * @param     y_type output type, F2C FAST-native derived type
 * @param     ierr error code
 * @param     map_msg error string
 * @see       {@link map_update_states()}
 */
MAP_EXTERNCALL void map_calc_output(float t,
                                    MAP_InputType_t* u_type,
                                    MAP_ParameterType_t* p_type,
                                    MAP_ContinuousStateType_t* x_type,
                                    MAP_DiscreteStateType_t* xd_type,
                                    MAP_ConstraintStateType_t* z_type,
                                    MAP_OtherStateType_t* other_type,
                                    MAP_OutputType_t* y_type,
                                    MAP_ERROR_CODE* ierr,
                                    char* map_msg);


/**
 * @brief     Deallocates all memory. Must be called after {@link map_init()}. This is called once. 
 * @param     u_type input type, F2C FAST-native derived type. 
 * @param     p_type parameter type, F2C FAST-native derived type
 * @param     x_type continuous-state type, F2C FAST-native derived type
 * @param     none void type. It exists because MAP does not use discrete types (it would fill this place otherwise)
 * @param     z_type constraint-state type, F2C FAST-native derived tpe
 * @param     other_type other-state type, F2C FAST-native derived type. This is core structure that houses all internal MAP states that cannot be mapped to Fortran
 * @param     y_type output type, F2C FAST-native derived type
 * @param     ierr error code
 * @param     map_msg error string
 * @see       {@link map_update_states()}
 */
MAP_EXTERNCALL void map_end(MAP_InputType_t* u_type,
                            MAP_ParameterType_t* p_type,
                            MAP_ContinuousStateType_t* x_type,
                            MAP_DiscreteStateType_t* xd_type,
                            MAP_ConstraintStateType_t* z_type,
                            MAP_OtherStateType_t* other_type,
                            MAP_OutputType_t* y_type,                                                                           
                            MAP_ERROR_CODE* ierr,
                            char* map_msg);

/**
 * @brief Set the name out the MAP summary output file. Does not need to be called; the default
 *        summary file name is 'outlist.map.sum'.
 * @param init_type initalization type, native C struct {@link InitializationData_t}
 * @param map_msg MAP error message
 * @param ierr MAP error code
 *
 * Example Fortran usage:
 * @code
 * ! Interface block declaration:
 * INTERFACE 
 *    SUBROUTINE mapextern_map_set_summary_file_name(fc_init, fc_char, fc_int) BIND(C,name='map_set_summary_file_name')
 *      IMPORT                                    
 *      IMPLICIT NONE                             
 *      TYPE(MAP_InitInputType_C) fc_init        
 *      CHARACTER(KIND=C_CHAR), DIMENSION(*) :: fc_char
 *      INTEGER(KIND=C_INT) :: fc_int            
 *    END SUBROUTINE mapextern_map_set_summary_file_name
 * END INTERFACE
 *
 *   ...
 *
 * ! access the function using this subroutine call: 
 * CALL mapextern_map_set_summary_file_name(InitInp%C_obj, ErrMsg, ErrStat)
 * @endcode
 * @todo: need to free summary_file_name. This is done in delete_all_init_data(...), should be called in Fortran routines
 */
MAP_EXTERNCALL void map_set_summary_file_name(MAP_InitInputType_t* init_type, char* map_msg, MAP_ERROR_CODE* ierr); 


/**
 * @brief Obtains the variable name array corresponding to the outputs selected in the MAP input file. For example,
 *        str_array can be:
 *        <pre>
 *        X[2]     H[1]     X[6]     H[3]     X[10]    H[5]     X[14]    H[7]
 *        </pre>
 * @param n number of header blocks. Should be proportional to the number of itms being output to the FAST output file
 * @param str_array the string being output.         
 * @param other_type Fortran other state derived type
 *
 * Example Fortran usage:
 * @code
 * ! Interface block declaration:
 * INTERFACE      
 *    SUBROUTINE mapextern_map_get_header_string(fc_int, fc_string, fc_other) BIND(C,name='map_get_header_string')   
 *      IMPORT                                 
 *      IMPLICIT NONE                          
 *      INTEGER(KIND=C_INT) :: fc_int          
 *      TYPE( MAP_OtherStateType_C ) fc_other                                                 
 *      TYPE(C_PTR), DIMENSION(FC_int) :: fc_string
 *    END SUBROUTINE mapextern_map_get_header_string
 * END INTERFACE                                                 
 *
 *   ...
 *
 * ! access the function using this subroutine call: 
 * CALL mapextern_map_get_header_string(num_header_str, hdr_str_ptrs, other%C_obj)
 * @endcode
 * @todo this should raise and error when count!=n
 */
MAP_EXTERNCALL void map_get_header_string(int* n, char** str_array, MAP_OtherStateType_t* other_type);


/**
 * @brief Obtains the units of the outputs passed back to the calling program. str_array can be:
 *        <pre>
 *        [m]     [N]     [m]     [N]     [m]     [N]     [m]     [N]   
 *        </pre>
 * @param n number of header blocks. Should be proportional to the number of itms being output to the FAST output file
 * @param str_array the string being output.         
 * @param other_type Fortran other state derived type
 *
 * Example Fortran usage:
 * @code
 * ! Interface block declaration:
 * INTERFACE
 *    SUBROUTINE mapextern_map_get_unit_string(fc_int, fc_string, fc_other) BIND(C,name='map_get_unit_string')          
 *      IMPORT                                     
 *      IMPLICIT NONE                              
 *      INTEGER(KIND=C_INT) :: fc_int              
 *      TYPE(MAP_OtherStateType_C) fc_other                                                 
 *      TYPE(C_PTR), DIMENSION(FC_int) :: fc_string
 *    END SUBROUTINE mapextern_map_get_unit_string        
 * END INTERFACE                                                 
 *
 *   ...
 * ! access the function using this subroutine call: 
 * CALL mapextern_map_get_header_string(num_header_str, unit_str_ptrs, other%C_obj)
 * @endcode
 * @todo this should raise and error when count!=n
 */
MAP_EXTERNCALL void map_get_unit_string(int* n, char** str_array ,MAP_OtherStateType_t* other_type);


/**
 * @brief   Allocate InitializationData
 * @details Called by {@link  map_create_init_type} to allocate memory for the iinitialization
 *          data type. The reason why a layer is added to the initialization data is due to 
 *          Fortran interactions. It is straighforward to pass 1D character arrays between
 *          Fortran and C instead of 2D arrays. 2D arrays would make more sense since multiple 
 *          lines from the MAP input file can be packed in one step. {@link MAP_InitInputType_t}
 *          in responsible for the 1D arrays. which are passed from Fortran to C. MAP then takes
 *          the 1D aray and packs it into InitializationData. This is used to subsequently 
 *          initialize the model. Structure is free'd by calling {@link MAP_InitInput_Delete}.
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  instance of the packed initialization strings (different from the FAST-required derived types)  
 */
MAP_EXTERNCALL InitializationData* MAP_InitInput_Create(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_InitInputType_t and InitializationData
 * @details Called to allocate memory for the initialzation data for both the Fortran
 *          derived data and internal state data. Following sucessful allocation, 
 *          {@link initialize_init_type_to_null} and {@link initialize_init_data_to_null}
 *          are both called to nullify data. If not called, memory errors results. This should 
 *          the first function called when interacting with MAP. This is a necessary function for
 *          interaction with python and C based programs
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  initialization input type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_InitInputType_t* map_create_init_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate Domain
 * @details Called by {@link  map_create_other_type} to allocate memory for the internal 
 *          state (model) data type. 'Other States', as FAST calls them, are states not 
 *          fitting a regular role as a parameter, constraint, input, ect. Other states
 *          contain information on the line connectivity matrix, how reference to poperties
 *          for each line, buoyancy properties of the nodes, ect. Deallocated using
 *          interaction with python and C based programs. Structure is free'd by calling
 *          {@link MAP_OtherState_Delete}.
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see     map_create_other_type()
 * @return  instance of the interal model struct (different from the FAST-required derived types)  
 */
MAP_EXTERNCALL Domain* MAP_OtherState_Create(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_OtherStateType_t and Domain
 * @details Called to allocate memory for the other states for both the Fortran
 *          derived data and internal state data. This is a necessary function for
 *          interaction with python and C based programs. The 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  other state type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_OtherStateType_t* map_create_other_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_InitOutputType_t 
 * @details Called to allocate memory for the initialization output type. The only obligation of
 *          this struct is to store the program version, necessary for FAST. This function is a
 *          necessary call for C and python, but can be ignored for Fortran if the MAP template 
 *          is followed (that is, ISO C Binding is followed in the mapping of Fortran type and C 
 *          structures). 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  initialization output type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_InitOutputType_t* map_create_initout_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_InputType_t
 * @details Called to allocate memory for the input type. The program inputs are the fairlead 
 *          displacement due to the motion of the vessel the cable are attached to. This function 
 *          is a necessary call for C and python, but can be ignored for Fortran if the MAP 
 *          template is followed (that is, ISO C Binding is followed in the mapping of Fortran 
 *          type and C structures). 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  input type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_InputType_t* map_create_input_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_ParameterType_t
 * @details Called to allocate memory for the parameter type. Parameters are time-invariant 
 *          constants, such as gravitational constant.  This function 
 *          is a necessary call for C and python, but can be ignored for Fortran if the MAP 
 *          template is followed (that is, ISO C Binding is followed in the mapping of Fortran 
 *          type and C structures). 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  parameter type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_ParameterType_t* map_create_parameter_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_ConstraintType_t
 * @details Called to allocate memory for the constraint type. Constraints are variables solved
 *          through an algebraic equation. This is fairlead end forces (H and V) and node positions.  
 *          This function is a necessary call for C and python, but can be ignored for Fortran if
 *          the MAP template is followed (that is, ISO C Binding is followed in the mapping of Fortran 
 *          type and C structures). 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  constraint type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_ConstraintStateType_t* map_create_constraint_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_OutputType_t
 * @details Called to allocate memory for the output type. IMPORTANT: this is different from the {@link OutList}.
 *          Output types are forces at the line fairlead only for lines connecting to the vessel. 
 *          This function is a necessary call for C and python, but can be ignored for Fortran if
 *          the MAP template is followed (that is, ISO C Binding is followed in the mapping of Fortran 
 *          type and C structures). 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  output type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_OutputType_t* map_create_output_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_ContinuousStateType_t
 * @details Called to allocate memory for the coninuous type. Not currently used, but it still is
 *          required to be allocated for FAST consistentcy. 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  continuous type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_ContinuousStateType_t* map_create_continuous_type(char* map_msg, MAP_ERROR_CODE* ierr);


#endif /* _MAPAPI_H */
