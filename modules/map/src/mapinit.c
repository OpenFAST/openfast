/*---------------------------------------------------------------
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
 ---------------------------------------------------------------*/


#include "mapinit.h"
#include "lineroutines.h"
#include "outputstream.h"


extern const char MAP_ERROR_STRING[][1024];


MAP_ERROR_CODE initialize_fortran_types(MAP_InputType_t* u_type, 
                                        MAP_ParameterType_t* p_type, 
                                        MAP_ContinuousStateType_t* x_type, 
                                        MAP_ConstraintStateType_t* z_type, 
                                        MAP_OtherStateType_t* other_type, 
                                        MAP_OutputType_t* y_type, 
                                        MAP_InitOutputType_t* initout_type)
{
  /* parameters are skipped for now; they are set in fortran since depth, 
   * gravity and sea density are set by glue code 
   */

  /* inputs */
  u_type->x = NULL;     u_type->x_Len = 0;
  u_type->y = NULL;     u_type->y_Len = 0;
  u_type->z = NULL;     u_type->z_Len = 0;

  /* continuous state */
  x_type->dummy=-999.9;

  /* constraint state */  
  z_type->H = NULL;     z_type->H_Len = 0;
  z_type->V = NULL;     z_type->V_Len = 0;
  z_type->x = NULL;     z_type->x_Len = 0;
  z_type->y = NULL;     z_type->y_Len = 0;
  z_type->z = NULL;     z_type->z_Len = 0;

  /* other state */
  other_type->H = NULL;     other_type->H_Len = 0;
  other_type->V = NULL;     other_type->V_Len = 0;
  other_type->Ha = NULL;    other_type->Ha_Len = 0;
  other_type->Va = NULL;    other_type->Va_Len = 0;
  other_type->x = NULL;     other_type->x_Len = 0;
  other_type->y = NULL;     other_type->y_Len = 0;
  other_type->z = NULL;     other_type->z_Len = 0;
  other_type->xa = NULL;    other_type->xa_Len = 0;
  other_type->ya = NULL;    other_type->ya_Len = 0;
  other_type->za = NULL;    other_type->za_Len = 0;
  other_type->Fx_connect = NULL;    other_type->Fx_connect_Len = 0;
  other_type->Fy_connect = NULL;    other_type->Fy_connect_Len = 0;
  other_type->Fz_connect = NULL;    other_type->Fz_connect_Len = 0;
  other_type->Fx_anchor = NULL;    other_type->Fx_anchor_Len = 0;
  other_type->Fy_anchor = NULL;    other_type->Fy_anchor_Len = 0;
  other_type->Fz_anchor = NULL;    other_type->Fz_anchor_Len = 0;

  /* outputs */
  y_type->Fx = NULL;              y_type->Fx_Len = 0;
  y_type->Fy = NULL;              y_type->Fy_Len = 0;
  y_type->Fz = NULL;              y_type->Fz_Len = 0;
  y_type->wrtOutput = NULL;       y_type->wrtOutput_Len = 0;
  y_type->WriteOutput = NULL;     y_type->WriteOutput_Len = 0;
  
  /* init outputs */
  initout_type->progName[0] = '\0' ;
  initout_type->version[0] = '\0';
  initout_type->compilingData[0] = '\0';
  initout_type->writeOutputHdr = NULL;     initout_type->writeOutputHdr_Len = 0;
  initout_type->writeOutputUnt = NULL;     initout_type->writeOutputUnt_Len = 0;

  return MAP_SAFE;
};


void initialize_init_type_to_null(MAP_InitInputType_t* init_type)
{
  /* initialize the native Fortran/C types */
  init_type->gravity = -999.9;
  init_type->sea_density = -999.9;
  init_type->depth = -999.9;
  init_type->file_name[0] = '\0';
  init_type->summary_file_name[0] = '\0';
  init_type->library_input_str[0] = '\0';
  init_type->node_input_str[0] = '\0';
  init_type->line_input_str[0] = '\0';
  init_type->option_input_str[0] = '\0';
};



void initialize_init_data_to_null(InitializationData* init_data)
{
  /* initialize the MAP initialization internal data strcture */
  init_data->library_input_string = bstrListCreate();
  init_data->node_input_string = bstrListCreate();
  init_data->line_input_string = bstrListCreate();
  init_data->solver_options_string = bstrListCreate();
  init_data->expanded_node_input_string = bstrListCreate();
  init_data->expanded_line_input_string = bstrListCreate();
};


void initialize_domain_to_null(Domain* domain)
{
  domain->MAP_SOLVE_TYPE = -999;
  domain->y_list = NULL; 
  initialize_inner_solve_data_defaults(&domain->inner_loop);    
  initialize_outer_solve_data_defaults(&domain->outer_loop);    
  initialize_vessel_to_null(&domain->vessel);    
  initialize_model_option_defaults(&domain->model_options);
};


size_t cable_line_meter(const void *el) 
{
  return sizeof(Line);
};


size_t cable_library_meter(const void *el) 
{
  return sizeof(CableLibrary);
};


size_t node_meter(const void *el) 
{
  return sizeof(Node);
};


size_t u_list_meter(const void *el) 
{
  return sizeof(ReferencePoint);
};


size_t vartype_meter(const void* el) 
{
  /* every line has the constant size of a rectangle structure */
  return sizeof(VarType);
};


size_t vartype_ptr_meter(const void* el) 
{
  /* every line has the constant size of a rectangle structure */
  return sizeof(VarTypePtr);
};


MAP_ERROR_CODE allocate_outlist(Domain* data, char* map_msg, MAP_ERROR_CODE* ierr)
{ 
  data->y_list = malloc(sizeof(OutputList)); 
  if (data->y_list==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_46);    
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


void initialize_inner_solve_data_defaults(InnerSolveAttributes* inner)
{
  inner->f_tol = 1e-6;
  inner->g_tol = 1e-6;
  inner->x_tol = 1e-6;
  inner->max_its = 500;
  inner->m = 2;
  inner->n = 2;
  inner->factor = 1.0E2;             
  inner->factor = 1.0E2;             
  inner->ldfjac = 2; 
  inner->mode = 1;             
  inner->nprint = 2;           
  inner->info = 0;             
};


void initialize_model_option_defaults(DomainOptions* options) 
{
  options->repeat_angle_size = 0;
  options->repeat_angle = NULL;
  options->integration_dt = 0.01;
  options->kb_lm = 3.0E6;
  options->cb_lm = 3.0E5;
  options->wave_kinematics = false;
  options->lm_model = false;
}; 


/*


*/

/* deallocated in free_outer_solve_data() */
void initialize_outer_solve_data_defaults(OuterSolveAttributes* outer) 
{
  outer->fd = BACKWARD_DIFFERENCE;
  outer->pg = false;
  outer->krylov_accelerator = false;
  outer->powell = false;
  outer->tol = 1e-6;
  outer->epsilon = 1e-3;
  outer->max_its = 500;
  outer->jac = NULL;
  outer->x = NULL;
  outer->b = NULL;
  outer->l = NULL;
  outer->u = NULL;
  outer->y = NULL;

  outer->max_krylov_its = 3;
  outer->AV = NULL;
  outer->V = NULL;
  outer->av = NULL;
  outer->C = NULL;
  outer->q = NULL;
  outer->w = NULL;
};


void initialize_vessel_to_null(Vessel* floater)
{
  floater->xi = NULL;
  floater->yi = NULL;
  floater->zi = NULL;
    
  floater->displacement.x.name = NULL;
  floater->displacement.y.name = NULL;
  floater->displacement.z.name = NULL;
  floater->ref_origin.x.name = NULL;
  floater->ref_origin.y.name = NULL;
  floater->ref_origin.z.name = NULL;
  floater->line_sum_force.fx.name = NULL;
  floater->line_sum_force.fy.name = NULL;
  floater->line_sum_force.fz.name = NULL;
  floater->orientation.phi.name = NULL;
  floater->orientation.the.name = NULL;
  floater->orientation.psi.name = NULL;
    
  floater->displacement.x.units = NULL;
  floater->displacement.y.units = NULL;
  floater->displacement.z.units = NULL;
  floater->ref_origin.x.units = NULL;
  floater->ref_origin.y.units = NULL;
  floater->ref_origin.z.units = NULL;
  floater->line_sum_force.fx.units = NULL;
  floater->line_sum_force.fy.units = NULL;
  floater->line_sum_force.fz.units = NULL;
  floater->orientation.phi.units = NULL;
  floater->orientation.the.units = NULL;
  floater->orientation.psi.units = NULL;

  floater->ref_origin.x.value = 0.0;
  floater->ref_origin.y.value = 0.0;
  floater->ref_origin.z.value = 0.0;
};


MAP_ERROR_CODE set_vessel(Vessel* floater, const MAP_InputType_t* u_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  int n = u_type->x_Len;

  MAP_BEGIN_ERROR_LOG;

  /* vessel displacement */
  success = set_vartype_float("[m]", "Vessel_X", -999, &floater->displacement.x, 0.0); CHECKERRQ(MAP_FATAL_68);
  success = set_vartype_float("[m]", "Vessel_Y", -999, &floater->displacement.y, 0.0); CHECKERRQ(MAP_FATAL_68);
  success = set_vartype_float("[m]", "Vessel_Z", -999, &floater->displacement.z, 0.0); CHECKERRQ(MAP_FATAL_68);
     
  /* vessel reference origin. When ==[0.0, 0.0, 0.0], then the reference origin is aligned with the SWL 
   * Note: this is commented because it over rides the run-time option 'REF_POSITION'. Instead, the ref position is 
   * initialized to zero in function void initialize_vessel_to_null(Vessel* floater)
   */
  //success = set_vartype_float("[m]", "Vessel_Xref", -999, &floater->ref_origin.x, 0.0); CHECKERRQ(MAP_FATAL_68);
  //success = set_vartype_float("[m]", "Vessel_Yref", -999, &floater->ref_origin.y, 0.0); CHECKERRQ(MAP_FATAL_68);
  //success = set_vartype_float("[m]", "Vessel_Zref", -999, &floater->ref_origin.z, 0.0); CHECKERRQ(MAP_FATAL_68);
    
  /* sum force of all fairleads connecte to the vessel */
  success = set_vartype_float("[N]", "Vessel_fx", -999, &floater->line_sum_force.fx, 0.0); CHECKERRQ(MAP_FATAL_68);
  success = set_vartype_float("[N]", "Vessel_fy", -999, &floater->line_sum_force.fy, 0.0); CHECKERRQ(MAP_FATAL_68);
  success = set_vartype_float("[N]", "Vessel_fz", -999, &floater->line_sum_force.fz, 0.0); CHECKERRQ(MAP_FATAL_68);
    
  /* orientation of the vessel. This is used as input from the user */
  success = set_vartype_float("[deg]", "Vessel_phi", -999, &floater->orientation.phi, 0.0); CHECKERRQ(MAP_FATAL_68);
  success = set_vartype_float("[deg]", "Vessel_the", -999, &floater->orientation.the, 0.0); CHECKERRQ(MAP_FATAL_68);
  success = set_vartype_float("[deg]", "Vessel_psi", -999, &floater->orientation.psi, 0.0); CHECKERRQ(MAP_FATAL_68);

  MAP_END_ERROR_LOG;

  floater->xi = malloc(n*sizeof(double));  
  floater->yi = malloc(n*sizeof(double));  
  floater->zi = malloc(n*sizeof(double));  

  if (floater->xi==NULL || floater->yi==NULL || floater->zi==NULL) {
    return MAP_FATAL;
  };
  
  for (i=0 ; i<n ; i++) {
    floater->xi[i] = u_type->x[i];
    floater->yi[i] = u_type->y[i];
    floater->zi[i] = u_type->z[i];
  };    
  return MAP_SAFE;
};


MAP_ERROR_CODE first_solve(Domain* domain, MAP_ParameterType_t* p_type, MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  
  if (domain->MAP_SOLVE_TYPE==MONOLITHIC) {
    success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr); /* @todo CHECKERRQ() */
  } else {
    success = node_solve_sequence(domain, p_type, u_type, z_type, other_type, 0.0, map_msg, ierr); /* @todo CHECKERRQ() */
  };
  
  MAP_RETURN_STATUS(success);
};


MAP_ERROR_CODE allocate_outer_solve_data(OuterSolveAttributes* ns, const int size, char* map_msg, MAP_ERROR_CODE* ierr)
{
  // int ret = 0;
  int i = 0;
  int j = 0;
  const int THREE = 3;  
  const int N = ns->max_krylov_its + 1;
  const int SIZE = THREE*size;

  ns->jac = malloc(SIZE*sizeof(double*));
  ns->l = malloc(SIZE*sizeof(double*));  
  ns->u = malloc(SIZE*sizeof(double*));  
  ns->x = malloc(SIZE*sizeof(double));
  ns->b = malloc(SIZE*sizeof(double));
  ns->y = malloc(SIZE*sizeof(double));  
  
  if (ns->jac==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_8);        
    return MAP_FATAL;
  };

  if (ns->x==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_8);        
    return MAP_FATAL;
  };

  if (ns->b==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_8);        
    return MAP_FATAL;
  };

  if (ns->l==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_8);        
    return MAP_FATAL;
  };

  if (ns->u==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_8);        
    return MAP_FATAL;
  };

  if (ns->y==NULL) {
    set_universal_error(map_msg, ierr, MAP_FATAL_8);        
    return MAP_FATAL;
  };

  /* final initialization to -999.9 */
  for(i=0 ; i<SIZE ; i++) {
    ns->jac[i] = malloc(SIZE*sizeof(double));    
    ns->l[i] = malloc(SIZE*sizeof(double));    
    ns->u[i] = malloc(SIZE*sizeof(double));    
    ns->x[i] = -999.9;
    ns->b[i] = -999.9;
    ns->y[i] = -999.9;

    if (ns->jac[i]==NULL) {
      set_universal_error(map_msg, ierr, MAP_FATAL_8);        
      return MAP_FATAL;
    };
    if (ns->l[i]==NULL) {
      set_universal_error(map_msg, ierr, MAP_FATAL_8);        
      return MAP_FATAL;
    };
    if (ns->u[i]==NULL) {
      set_universal_error(map_msg, ierr, MAP_FATAL_8);        
      return MAP_FATAL;
    };
    for(j=0 ; j<SIZE ; j++) {
      ns->jac[i][j] = -999.9;
      ns->l[i][j] = -999.9;   
      ns->u[i][j] = -999.9;
    };
  };
  
  if (ns->krylov_accelerator) { /* only allocated if  Krylov accelerator algorimth is invoked */    
    ns->AV = malloc(SIZE*sizeof(double*));
    ns->V = malloc(SIZE*sizeof(double*));
    ns->av = malloc(SIZE*N*sizeof(double));
    ns->C = malloc(SIZE*sizeof(double));  
    ns->q = malloc(SIZE*sizeof(double));  
    ns->w = malloc(SIZE*sizeof(double));  

    if (ns->AV==NULL) {
      set_universal_error(map_msg, ierr, MAP_FATAL_8);        
      return MAP_FATAL;
    };

    if (ns->V==NULL) {
      set_universal_error(map_msg, ierr, MAP_FATAL_8);        
      return MAP_FATAL;
    };

    if (ns->C==NULL) {
      set_universal_error(map_msg, ierr, MAP_FATAL_8);        
      return MAP_FATAL;
    };

    if (ns->q==NULL) {
      set_universal_error(map_msg, ierr, MAP_FATAL_8);        
      return MAP_FATAL;
    };

    if (ns->w==NULL) {
      set_universal_error(map_msg, ierr, MAP_FATAL_8);        
      return MAP_FATAL;
    };

    /* final initialization to -999.9 */
    for(i=0 ; i<SIZE ; i++) {
      ns->AV[i] = malloc(N*sizeof(double));    
      ns->V[i] = malloc(N*sizeof(double));    
      if (ns->AV[i]==NULL) {
        set_universal_error(map_msg, ierr, MAP_FATAL_8);        
        return MAP_FATAL;
      };
      if (ns->V[i]==NULL) {
        set_universal_error(map_msg, ierr, MAP_FATAL_8);        
        return MAP_FATAL;
      };
      ns->C[i] = -999.9;
      ns->q[i] = -999.9;
      ns->w[i] = -999.9;
      for(j=0 ; j<N ; j++) {
        ns->AV[i][j] = -999.9;
        ns->V[i][j] = -999.9;
      };
    };
  };

  return MAP_SAFE;
};


MAP_ERROR_CODE check_help_flag(bstring list)
{
  MAP_ERROR_CODE success = 0;

  success = biseqcstrcaseless(list,"HELP"); /* string compare */
  if (success) { 
    print_help_to_screen();
  }; 
  return MAP_SAFE;
};


MAP_ERROR_CODE check_inner_f_tol_flag(struct bstrList* list, double* ftol)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"INNER_FTOL"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = (char*)list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *ftol = (double)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_inner_g_tol_flag(struct bstrList* list, double* gtol)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"INNER_GTOL"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = (char*)list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *gtol = (double)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_inner_x_tol_flag(struct bstrList* list, double* xtol)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"INNER_XTOL"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = (char*)list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *xtol = (double)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};



MAP_ERROR_CODE check_inner_max_its_flag(struct bstrList* list, int* max_its)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"INNER_MAX_ITS"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = (char*)list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *max_its = (int)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_outer_max_its_flag(struct bstrList* list, int* max_its)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"OUTER_MAX_ITS"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = (char*)list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *max_its = (int)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_outer_tol_flag(struct bstrList* list, double* outer_tol)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"OUTER_TOL"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = (char*)list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *outer_tol = (double)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_outer_epsilon_flag(struct bstrList* list, double* epsilon)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"OUTER_EPSILON"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = (char*)list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *epsilon = (double)atof(word);
          return MAP_SAFE;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_integration_dt_flag(struct bstrList* list, double* dt)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"INTEGRATION_DT"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = (char*)list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *dt = (double)atof(word);
          return MAP_WARNING;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_kb_default_flag(struct bstrList* list, double* kb)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"KB_DEFAULT"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = (char*)list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *kb = (double)atof(word);
          return MAP_WARNING;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_krylov_accelerator_flag(struct bstrList* list, OuterSolveAttributes* solver)
{
  int n = 0;
  int success = 0;
  // int next = 0; 
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"KRYLOV_ACCELERATOR"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) {
    solver->krylov_accelerator = true;
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = (char*)list->entry[n+1]->data;
        if (is_numeric(word)) {         
          solver->max_krylov_its = (int)atoi(word);          
          word = NULL;
          return MAP_SAFE;        
        };
      } else { /* no trailing integer in the MAP input file */
        word = NULL;
        return MAP_WARNING;
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_powell_flag(struct bstrList* list, OuterSolveAttributes* solver)
{
  int success = 0;

  success = biseqcstrcaseless(list->entry[0],"POWELL"); /* string compare */
  if (success) {
    solver->powell = true;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE check_cb_default_flag(struct bstrList* list, double* cb)
{
  int success = 0;
  int n = 0;
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"CB_DEFAULT"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) { 
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = (char*)list->entry[n+1]->data;
        if (is_numeric(word)) { 
          *cb = (double)atof(word);
          return MAP_WARNING;
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_outer_bd_flag(struct bstrList* list, FdType* bd)
{
  int success = 0;

  success = biseqcstrcaseless(list->entry[0],"OUTER_BD"); /* string compare */
  if (success) {
    *bd = BACKWARD_DIFFERENCE;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE check_outer_cd_flag(struct bstrList* list, FdType* cd)
{
  int success = 0;

  success = biseqcstrcaseless(list->entry[0],"OUTER_CD"); /* string compare */
  if (success) {
    *cd = CENTRAL_DIFFERENCE;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE check_outer_fd_flag(struct bstrList* list, FdType* fd)
{
  int success = 0;

  success = biseqcstrcaseless(list->entry[0],"OUTER_FD"); /* string compare */
  if (success) {
    *fd = FORWARD_DIFFERENCE;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE check_wave_kinematics_flag(struct bstrList* list, bool* wave)
{
  int success = 0;

  success = biseqcstrcaseless(list->entry[0],"WAVE_KINEMATICS"); /* string compare */
  if (success) {
    *wave = false;
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE check_lm_model_flag(struct bstrList* list, bool* lm)
{
  int success = 0;

  success = biseqcstrcaseless(list->entry[0],"LM_MODEL"); /* string compare */
  if (success) {
    *lm = true;
    return MAP_SAFE;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE check_pg_cooked_flag(struct bstrList* list, OuterSolveAttributes* solver)
{
  int n = 0;
  int success = 0;
  int next = 0; 
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"PG_COOKED"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) {
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = (char*)list->entry[n+1]->data;
        if (is_numeric(word)) {         
          if (!next) {
            solver->d = (double)atof(word);
            next++;
          } else {
            solver->ds = (double)atof(word);
            solver->pg = true;
            return MAP_SAFE;
          };
        } else {
          return MAP_FATAL;
        };
      };
      n++;
    };
  }; 
  word = NULL;
  if (!next) {
    return MAP_SAFE;
  } else {
    return MAP_WARNING;
  };
};


MAP_ERROR_CODE check_repeat_flag(struct bstrList* list, DomainOptions* options)
{
  double* more_angles = NULL;
  char* current = NULL;
  int success = 0;
  int n = 0; /* word interator in the list */
  int i = 0; /* synonym for repeat_angle_size */

  success = biseqcstrcaseless(list->entry[0],"REPEAT"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) {
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        current = (char*)list->entry[n+1]->data;
        i = options->repeat_angle_size;
        more_angles = realloc(options->repeat_angle, (i+1)*sizeof(double));
        if (more_angles) {
          options->repeat_angle = more_angles;
          if (is_numeric(current)) { /* add repeat angle if word is numeric */
            options->repeat_angle[i] = atof(current);
            options->repeat_angle_size++;
          } else { 
            MAPFREE(more_angles);
            return MAP_FATAL;
          };
        } else {
          MAPFREE(more_angles);
          return MAP_FATAL;
        };        
      };
      n++;
    };
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE check_ref_position_flag(struct bstrList* list, Point* ref_position)
{
  int n = 0;
  int success = 0;
  int next = 0; 
  const char* word = NULL;

  success = biseqcstrcaseless(list->entry[0],"REF_POSITION"); /* string compare */
  if (success==BSTR_ERR) {
    return MAP_FATAL;
  } else if (success) {
    while (n<list->qty-1) { /* iterating through all strings */      
      if (list->entry[n+1]->slen) { /* if the string length is not 0 */
        word = (char*)list->entry[n+1]->data;
        if (is_numeric(word)) {         
          if (!next) {
            ref_position->x.value = (double)atof(word);
            next++;
          } else if (next==1) {
            ref_position->y.value = (double)atof(word);
            next++;
          } else {
            ref_position->z.value = (double)atof(word);
            return MAP_SAFE;
          };
        };
      };
      n++;
    };
    return MAP_WARNING;
  };
  word = NULL;
  return MAP_SAFE;
};


MAP_ERROR_CODE check_uncaught_flag(struct bstrList* list)
{
  if (biseqcstrcaseless(list->entry[0],"")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"HELP")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"INNER_FTOL")) { 
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"INNER_GTOL")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"INNER_XTOL")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"INNER_MAX_ITS")) { 
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"OUTER_MAX_ITS")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"OUTER_TOL")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"OUTER_EPSILON")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"INTEGRATION_DT")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"KB_DEFAULT")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"CB_DEFAULT")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"OUTER_CD")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"OUTER_BD")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"OUTER_FD")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"WAVE_KINEMATICS")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"LM_MODEL")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"POWELL")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"PG_COOKED")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"KRYLOV_ACCELERATOR")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"REPEAT")) {
    return MAP_SAFE;
  } else if (biseqcstrcaseless(list->entry[0],"REF_POSITION")) {
    return MAP_SAFE;
  }; 
  return MAP_WARNING;
};


MAP_ERROR_CODE set_library_diameter(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric((char*)word->data)) { 
    library_ptr->diam = (double)atof((char*)word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_mass_density(bstring word, CableLibrary* library_ptr)
{  
  if (is_numeric((char*)word->data)) { 
    library_ptr->mass_density = (double)atof((char*)word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_ea(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric((char*)word->data)) { 
    library_ptr->EA = (double)atof((char*)word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_cb(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric((char*)word->data)) { 
    library_ptr->cb = (double)atof((char*)word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_internal_damping(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric((char*)word->data)) { 
    library_ptr->cd_i = (double)atof((char*)word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_added_mass_coefficient(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric((char*)word->data)) { 
    library_ptr->ca = (double)atof((char*)word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_cross_flow_drag_coefficient(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric((char*)word->data)) { 
    library_ptr->cd_n = (double)atof((char*)word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE set_library_tangent_drag_coefficient(bstring word, CableLibrary* library_ptr)
{
  if (is_numeric((char*)word->data)) { 
    library_ptr->cd_t = (double)atof((char*)word->data);
  } else { 
    return MAP_FATAL;
  };  
  return MAP_SAFE;
};



MAP_ERROR_CODE set_model_options_list(Domain* domain, InitializationData* init_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  const int n_lines = (init_data->solver_options_string->qty)-1;
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 

  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   
  for (i=0 ; i<=n_lines ; i++) { 
    parsed = bsplits(init_data->solver_options_string->entry[i], &tokens);

    MAP_BEGIN_ERROR_LOG;

    success = check_help_flag(parsed->entry[0]); CHECKERRQ(MAP_FATAL_85);
    success = check_inner_f_tol_flag(parsed, &domain->inner_loop.f_tol); CHECKERRK(MAP_ERROR_2);
    success = check_outer_max_its_flag(parsed, &domain->outer_loop.max_its); CHECKERRK(MAP_ERROR_3);
    success = check_inner_max_its_flag(parsed, &domain->inner_loop.max_its); CHECKERRK(MAP_ERROR_4);
    success = check_inner_g_tol_flag(parsed, &domain->inner_loop.g_tol); CHECKERRK(MAP_ERROR_9);
    success = check_inner_x_tol_flag(parsed, &domain->inner_loop.x_tol); CHECKERRK(MAP_ERROR_10);
    success = check_outer_tol_flag(parsed, &domain->outer_loop.tol); CHECKERRK(MAP_ERROR_3);
    success = check_outer_epsilon_flag(parsed, &domain->outer_loop.epsilon); CHECKERRK(MAP_ERROR_3);
    success = check_integration_dt_flag(parsed, &domain->model_options.integration_dt); CHECKERRK(MAP_ERROR_15); 
    success = check_kb_default_flag(parsed, &domain->model_options.kb_lm); CHECKERRK(MAP_ERROR_16); 
    success = check_cb_default_flag(parsed, &domain->model_options.cb_lm); CHECKERRK(MAP_ERROR_17); 
    success = check_outer_bd_flag(parsed, &domain->outer_loop.fd);
    success = check_outer_cd_flag(parsed, &domain->outer_loop.fd);
    success = check_outer_fd_flag(parsed, &domain->outer_loop.fd);      
    success = check_wave_kinematics_flag(parsed, &domain->model_options.wave_kinematics); CHECKERRK(MAP_WARNING_10);
    success = check_lm_model_flag(parsed, &domain->model_options.lm_model); CHECKERRK(MAP_WARNING_11);
    success = check_pg_cooked_flag(parsed, &domain->outer_loop); CHECKERRK(MAP_WARNING_8);
    success = check_krylov_accelerator_flag(parsed, &domain->outer_loop); CHECKERRK(MAP_WARNING_12);
    success = check_powell_flag(parsed, &domain->outer_loop); CHECKERRK(MAP_WARNING_14);
    success = check_repeat_flag(parsed, &domain->model_options); CHECKERRQ(MAP_FATAL_34);
    success = check_ref_position_flag(parsed, &domain->vessel.ref_origin); CHECKERRQ(MAP_FATAL_36);
    success = check_uncaught_flag(parsed);       
    if (success) {
      set_universal_error_with_message(map_msg, ierr, MAP_WARNING_1, "word: <%s>", parsed->entry[0]->data);
    };

    MAP_END_ERROR_LOG;
   
    success = bstrListDestroy(parsed);
  };

  /* throw error if PG_COOKED and KRYLOV_ACCELERATOR are simultaneously set in the input file */
  if (domain->outer_loop.pg && domain->outer_loop.krylov_accelerator) {
    set_universal_error(map_msg, ierr, MAP_WARNING_13);
    domain->outer_loop.krylov_accelerator = false;
  };

  /* throw error if PG_COOKED and KRYLOV_ACCELERATOR are simultaneously set in the input file */
  if ( (domain->outer_loop.powell && domain->outer_loop.krylov_accelerator) ||
       (domain->outer_loop.powell && domain->outer_loop.pg)) {
    set_universal_error(map_msg, ierr, MAP_WARNING_14);
    domain->outer_loop.krylov_accelerator = false;
  };

  MAP_RETURN_STATUS(*ierr);
};


MAP_ERROR_CODE reset_cable_library(CableLibrary* library_ptr)
{
  library_ptr->diam = 0.0;
  library_ptr->mass_density = 0.0;
  library_ptr->EA = 0.0;          
  library_ptr->omega = 0.0;       
  library_ptr->a = 0.0;           
  library_ptr->cb = 0.0;          
  library_ptr->cd_i = 0.0;    
  library_ptr->ca = 0.0;      
  library_ptr->cd_n = 0.0; 
  library_ptr->cd_t = 0.0; 
  // if (library_ptr->label) {
  //   success = bdestroy(library_ptr->label);
  // };
  return MAP_SAFE;
};


MAP_ERROR_CODE set_cable_library_list(Domain* domain, InitializationData* init_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  int n = 0;
  int next = 0; 
  const int n_lines = (init_data->library_input_string->qty)-1;
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  CableLibrary new_cable_library;
  CableLibrary* library_iter = NULL;
  
  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   
  success = reset_cable_library(&new_cable_library);

  for (i=0 ; i<=n_lines ; i++) { 
    list_append(&domain->library, &new_cable_library);
    library_iter = (CableLibrary*)list_get_at(&domain->library, i);

    parsed = bsplits(init_data->library_input_string->entry[i], &tokens);
    n = 0;
    next = 0;

    MAP_BEGIN_ERROR_LOG;  

    while (n<parsed->qty-1) { /* iterating through all strings */              
      if (parsed->entry[n]->slen) { /* if the string length is not 0 */
        if (next==0) {
          library_iter->label = bstrcpy(parsed->entry[n]);                         
          next++;
        } else if (next==1) {
          success = set_library_diameter(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_12);
          next++;            
        } else if (next==2) {
          success = set_library_mass_density(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_13);
          next++;
        } else if (next==3) {
          success = set_library_ea(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_14);
          next++;
        } else if (next==4) {
          success = set_library_cb(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_15);
          next++;
        } else if (next==5) {
          success = set_library_internal_damping(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_81);
          next++;
        } else if (next==6) {
          success = set_library_added_mass_coefficient(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_82);
          next++;
        } else if (next==7) {
          success = set_library_cross_flow_drag_coefficient(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_83);
          next++;
        } else if (next==8) {
          success = set_library_tangent_drag_coefficient(parsed->entry[n], library_iter); CHECKERRQ(MAP_FATAL_84);
          next++;
        };
      };
      n++;
    };

    MAP_END_ERROR_LOG;   

    success = bstrListDestroy(parsed);
  };
  MAP_RETURN_STATUS(*ierr);
};


MAP_ERROR_CODE initialize_cable_library_variables(Domain* domain, MAP_ParameterType_t* p_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  double radius = 0.0;
  double area = 0.0;
  double mu = 0.0;
  double rho_fluid = 0.0; 
  const double g = p_type->g;
  const double PI = 3.14159264;
  CableLibrary* library_iter = NULL;

  list_iterator_start(&domain->library); /* starting an iteration "session" */
  while ( list_iterator_hasnext(&domain->library)) { /* tell whether more values available */ 
    library_iter = (CableLibrary*)list_iterator_next(&domain->library);
    radius = library_iter->diam/2;
    area = PI*pow(radius,2);
    mu = library_iter->mass_density;
    rho_fluid = p_type->rho_sea;
    library_iter->omega = g*(mu-area*rho_fluid);

    library_iter->a = area;
    if (fabs(library_iter->omega)<=1.0) {
      set_universal_error_with_message(map_msg, ierr, MAP_WARNING_5, 
                                       "omega = %f <= 1.0", library_iter->omega);
    };
  };
  list_iterator_stop(&domain->library); /* ending the iteration "session" */    
  
  if (fabs(library_iter->omega)<=1.0E-3) {
    return MAP_FATAL;
  }
  /* end read */
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_number(const int n_line, bstring line)
{
  bstring current_entry = NULL;
  int ret = 0;

  current_entry = bformat("%d   ", n_line);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);               
  return MAP_SAFE;
};

MAP_ERROR_CODE expand_node_type(const char* word, bstring line)
{
  bstring current_entry = NULL;
  int ret = 0;

  current_entry = bformat("%s   ", word);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_position_x(double* x, const char* word)
{
  if (is_numeric(word)) { /* if number is numeric */
    *x = (double)atof(word);                
  } else if (word[0]=='#') {    
    *x = (double)atof(remove_first_character(word));
  } else {
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_position_y(double* y, const char* word)
{
  if (is_numeric(word)) { /* if number is numeric */
    *y = (double)atof(word);                
  } else if (word[0]=='#') {    
    *y = (double)atof(remove_first_character(word));    
  } else {
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_position_z(Vector* position, const double angle, const double x, const double y, const char* word, bstring line)
{
  bstring current_entry = NULL;
  bstring compared_word = NULL;
  int ret = 0;

  position->x =  x*cos(angle) + y*sin(angle);
  position->y = -x*sin(angle) + y*cos(angle);                    
  compared_word = bformat("%s", word);
  if (is_numeric(word)) { /* if number is numeric */
    position->z = (double)atof(word);                
    current_entry = bformat("%1.4f   %1.4f   %1.4f   ",position->x, position->y, position->z);
  } else if (word[0]=='#') {    
    position->z = (double)atof(remove_first_character(word));    
    current_entry = bformat("#%1.4f   #%1.4f   #%1.4f   ",position->x, position->y, position->z);
  } else if (biseqcstrcaseless(compared_word,"DEPTH")) {
    current_entry = bformat("%1.4f   %1.4f   depth   ",position->x, position->y, position->z);
  } else {
    bdestroy(compared_word);
    return MAP_FATAL;
  };

  bdestroy(compared_word);  
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_mass(const char* word, bstring line)
{
  bstring current_entry = NULL;
  int ret = 0;

  current_entry = bformat("%s   ", word);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_buoyancy(const char* word, bstring line)
{
  bstring current_entry = NULL;
  int ret = 0;

  current_entry = bformat("%s   ", word);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_force_x(double* fx, const char* word)
{
  if (is_numeric(word)) { /* if number is numeric */
    *fx = (double)atof(word);                
  } else if (word[0]=='#') { /* if the nuymber is iterated */
    if (is_numeric(remove_first_character(word))) { 
      *fx = (double)atof(remove_first_character(word));
    } else { /* in this case, it is presumed the force is just '#' */
      *fx = -999.9;                  
    };
  } else {
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_force_y(double* fy, const char* word)
{
  if (is_numeric(word)) { /* if number is numeric */
    *fy = (double)atof(word);                
  } else if (word[0]=='#') { /* if the nuymber is iterated */
    if (is_numeric(remove_first_character(word))) { 
      *fy = (double)atof(remove_first_character(word));
    } else { /* in this case, it is presumed the force is just '#' */
      *fy = -999.9;                  
    };
  } else {
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_node_force_z(Vector* force, const double angle, const double fx, const double fy, const char* word, bstring line)
{
  bstring current_entry = NULL;
  int ret = 0;

  force->x =  fx*cos(angle) + fy*sin(angle);
  force->y = -fx*sin(angle) + fy*cos(angle);                    

  if (is_numeric(word)) { /* if number is numeric */
    force->z = (double)atof(word);                
    current_entry = bformat("%1.4f   %1.4f   %1.4f\n",force->x, force->y, force->z);              
  } else if (word[0]=='#') { /* if the nuymber is iterated */
    if (is_numeric(remove_first_character(word))) { 
      force->z = (double)atof(remove_first_character(word));
      current_entry = bformat("#%1.4f   #%1.4f   #%1.4f\n",force->x, force->y, force->z);              
    } else { /* in this case, it is presumed the force is just '#' */
      force->z = -999.9;                  
      current_entry = bformat("#   #   #\n");
    };
  } else {
    ret = bconcat(line, current_entry);
    ret = bdestroy(current_entry);               
    return MAP_FATAL;
  };

  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);               
  return MAP_SAFE;
};


MAP_ERROR_CODE repeat_nodes(Domain* domain, InitializationData* init_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  int j = 0;
  int next = 0; 
  int i_parsed = 0;
  int n_line = 0;
  const int num_repeat = domain->model_options.repeat_angle_size; 
  const int num_node = init_data->node_input_string->qty;
  const int n = (num_node)*(num_repeat+1);
  const char* word = NULL;  
  double x_position = 0.0;
  double y_position = 0.0;
  double x_force = 0.0;
  double y_force = 0.0;
  double current_angle = 0.0;
  bstring line = bformat("");            
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  Vector force;     
  Vector position;
  Node new_node; 

  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   

  force.x = 0.0;
  force.y = 0.0;
  force.z = 0.0;
  position.x = 0.0;
  position.y = 0.0;
  position.z = 0.0;
   
  /* allocate space needed to expand the number of node lines */
  success = bstrListAlloc(init_data->expanded_node_input_string, n+1); 
  init_data->expanded_node_input_string->qty = 0;
  
  for(i=0 ; i<num_node ; i++) {     
    init_data->expanded_node_input_string->entry[i] = bfromcstr((char*)init_data->node_input_string->entry[i]->data);// bstrcpy(init_data->node_input_string->entry[i]);
    init_data->expanded_node_input_string->qty++;
  };

  for(i=0 ; i<num_repeat ; i++) { /* this is skipped if not repeat angles are declared */
    for(j=0 ; j<num_node ; j++) { 
      success = reset_node(&new_node);  
      n_line = (i+1)*num_node + j;
      current_angle = domain->model_options.repeat_angle[i]*(DEG2RAD);
      parsed = bsplits(init_data->node_input_string->entry[j], &tokens);
      next = 0;
      i_parsed = 0;

      MAP_BEGIN_ERROR_LOG;  

      while (i_parsed<parsed->qty-1) { /* iterating through all strings */              
        if (parsed->entry[i_parsed]->slen) { /* if the string length is not 0 */
          word = (char*)parsed->entry[i_parsed]->data;      
          if (next==0) {
            success = expand_node_number(n_line+1, line);/* @todo: checkerrq */
            next++;
          } else if (next==1) {
            success = expand_node_type(word, line);/* @todo: checkerrq */
            next++;
          } else if (next==2) {
            success = expand_node_position_x(&x_position, word);/* @todo: checkerrq */
            next++;
          } else if (next==3) {
            success = expand_node_position_y(&y_position, word);/* @todo: checkerrq */
            next++;
          } else if (next==4) {
            success = expand_node_position_z(&position, current_angle, x_position, y_position, word, line);/* @todo: checkerrq */
            next++;
          } else if (next==5) { /* node mass */
            success = expand_node_mass(word, line);/* @todo: checkerrq */
            next++;
          } else if (next==6) { /* node buoyancy */
            success = expand_node_buoyancy(word, line);/* @todo: checkerrq */
            next++;
          } else if (next==7) {
            success = expand_node_force_x(&x_force, word);/* @todo: checkerrq */
            next++;
          } else if (next==8) {
            success = expand_node_force_y(&y_force, word);/* @todo: checkerrq */
            next++;
          } else if (next==9) {
            success = expand_node_force_z(&force, current_angle, x_force, y_force, word, line);/* @todo: checkerrq */
            next++;
          };
        };
        i_parsed++;
      };
      init_data->expanded_node_input_string->qty++;
      init_data->expanded_node_input_string->entry[n_line] = bstrcpy(line);
      success = bassigncstr(line, "");

      MAP_END_ERROR_LOG;   

      success = bstrListDestroy(parsed);
    };  
  };
  success = bdestroy(line);               

  MAP_RETURN_STATUS(*ierr);
};


MAP_ERROR_CODE expand_line_number(const int n_line, bstring line)
{
  bstring current_entry = NULL;
  int ret = 0;
              
  current_entry = bformat("%d   ", n_line);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);               
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_line_property_name(const char* word, bstring line)
{
  int ret = 0;
  bstring current_entry = NULL;

  current_entry = bformat("%s   ", word);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_line_length(const char* word, bstring line)
{
  int ret = 0;
  bstring current_entry = NULL;
  
  current_entry = bformat("%s   ", word);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_line_anchor_number(const char* word, const int index, const int n, bstring line)
{
  int ret = 0;
  bstring current_entry = NULL;
  int anchor_num = 0;

  if (is_numeric(word)) {
    anchor_num = (int)atoi(word);
  } else {
    return MAP_SAFE;
  };  
  current_entry = bformat("%d   ", (index+1)*n+anchor_num);
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE expand_line_fairlead_number(const char* word, const int index, const int n, bstring line)
{
  int ret = 0;
  bstring current_entry = NULL;
  int fairlead_num = 0;

  if (is_numeric(word)) {
    fairlead_num = (int)atoi(word);
  } else {
    return MAP_SAFE;
  };  
  current_entry = bformat("%d   ", (index+1)*n+fairlead_num);
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};

MAP_ERROR_CODE expand_line_flag(const char* word, bstring line)
{
  int ret = 0;
  bstring current_entry = NULL;
  
  current_entry = bformat("%s   ", word);              
  ret = bconcat(line, current_entry);
  ret = bdestroy(current_entry);
  return MAP_SAFE;
};


MAP_ERROR_CODE repeat_lines(Domain* domain, InitializationData* init_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  int j = 0;
  int next = 0; 
  int i_parsed = 0;
  int n_line = 0;
  const int num_repeat = domain->model_options.repeat_angle_size; 
  const int num_line = init_data->line_input_string->qty;
  const int num_node = init_data->node_input_string->qty;
  const int n = (num_line)*(num_repeat+1);
  const char* word = NULL;  
  double current_angle = 0.0;
  bstring line = bformat("");            
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  Line new_line;

  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   

  /* allocate space needed to expand the number of node lines */
  success = bstrListAlloc(init_data->expanded_line_input_string, n+1); 
  init_data->expanded_line_input_string->qty = 0;
  
  for(i=0 ; i<num_line ; i++) {     
    init_data->expanded_line_input_string->entry[i] = bfromcstr((char*)init_data->line_input_string->entry[i]->data);
    init_data->expanded_line_input_string->qty++;
  };
  

  for(i=0 ; i<num_repeat ; i++) { /* this is skipped if not repeat angles are declared */
    for(j=0 ; j<num_line ; j++) { 
      success = reset_line(&new_line);  
      n_line = (i+1)*num_line + j;
      current_angle = domain->model_options.repeat_angle[i]*(DEG2RAD);
      parsed = bsplits(init_data->line_input_string->entry[j], &tokens);
      next = 0;
      i_parsed = 0;

      MAP_BEGIN_ERROR_LOG;
  
      while (i_parsed<parsed->qty-1) { /* iterating through all strings */              
        if (parsed->entry[i_parsed]->slen) { /* if the string length is not 0 */
          word = (char*)parsed->entry[i_parsed]->data;
          if (next==0) {
            success = expand_line_number(n_line+1, line);
            next++;
          } else if (next==1) {
            success = expand_line_property_name(word, line);
            next++;
          } else if (next==2) {
            success = expand_line_length(word, line);
            next++;
          } else if (next==3) {
            success = expand_line_anchor_number(word, i, num_node, line);
            next++;
          } else if (next==4) {
            success = expand_line_fairlead_number(word, i, num_node, line);
            next++;
          } else {
            success = expand_line_flag(word, line);
            next++;
          };
        };
        i_parsed++;
      };
      init_data->expanded_line_input_string->qty++;
      init_data->expanded_line_input_string->entry[n_line] = bstrcpy(line);
      success = bassigncstr(line, "");

      MAP_END_ERROR_LOG;      

      success = bstrListDestroy(parsed);
    };
  };
  success = bdestroy(line);               

  MAP_RETURN_STATUS(*ierr);
};



MAP_ERROR_CODE allocate_types_for_nodes(MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, Domain* domain, struct bstrList* node_input_string, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int i_parsed = 0;
  int next = 0;
  int fix_num = 0;
  int vessel_num = 0;
  int connect_num = 0;
  MAP_ERROR_CODE success = MAP_SAFE;
  const int num_nodes = node_input_string->qty;
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  
  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   

   /* First need to run loop to find the number of inputs, outputs, constraints.
    * Basically we are just counting the number of VESSEL, FIX, and CONNECT nodes     
    */
  for(i=0 ; i<num_nodes ; i++) {             
    i_parsed = 0;
    next = 0;
    parsed = bsplits(node_input_string->entry[i], &tokens);

    MAP_BEGIN_ERROR_LOG;  

    while (i_parsed<parsed->qty-1) { /* iterating through all strings */              
      if (parsed->entry[i_parsed]->slen) { /* if the string length is not 0 */
        if (next==1) {
          if (biseqcstrcaseless(parsed->entry[i_parsed],"FIX")) {
            fix_num++;
            break; /* break the while-loop because the agenda is reached */
          } else if (biseqcstrcaseless(parsed->entry[i_parsed],"CONNECT")) {
            connect_num++;
            break; /* break the while-loop because the agenda is reached */
          } else if (biseqcstrcaseless(parsed->entry[i_parsed],"VESSEL")) {
            vessel_num++;
            break; /* break the while-loop because the agenda is reached */
          } else {
            set_universal_error_with_message(map_msg, ierr, MAP_FATAL_25, "Value: <%s>", parsed->entry[i_parsed]->data);
          };
        };          
        next++;
      };
      i_parsed++;
    };

    MAP_END_ERROR_LOG;

    success = bstrListDestroy(parsed);
  }; 

  other_type->x_Len = fix_num;
  other_type->y_Len = fix_num;
  other_type->z_Len = fix_num;
  other_type->x = malloc(other_type->x_Len*sizeof(double));
  other_type->y = malloc(other_type->y_Len*sizeof(double));
  other_type->z = malloc(other_type->z_Len*sizeof(double));

  /* If the node is VESSEL, then the applied force is an output state. Otherwise, 
   * it has to be an other state because it can't be associated with any other type. 
   * This is what is done below.
   */
  other_type->Fx_connect_Len = connect_num;
  other_type->Fy_connect_Len = connect_num;
  other_type->Fz_connect_Len = connect_num;
  other_type->Fx_connect = malloc(other_type->Fx_connect_Len*sizeof(double));
  other_type->Fy_connect = malloc(other_type->Fy_connect_Len*sizeof(double));
  other_type->Fz_connect = malloc(other_type->Fz_connect_Len*sizeof(double));

  other_type->Fx_anchor_Len = fix_num;
  other_type->Fy_anchor_Len = fix_num;
  other_type->Fz_anchor_Len = fix_num;
  other_type->Fx_anchor = malloc(other_type->Fx_anchor_Len*sizeof(double));
  other_type->Fy_anchor = malloc(other_type->Fy_anchor_Len*sizeof(double));
  other_type->Fz_anchor = malloc(other_type->Fz_anchor_Len*sizeof(double));
  
  z_type->x_Len = connect_num;          
  z_type->y_Len = connect_num;          
  z_type->z_Len = connect_num;          
  z_type->x = malloc(z_type->x_Len*sizeof(double));
  z_type->y = malloc(z_type->y_Len*sizeof(double));
  z_type->z = malloc(z_type->z_Len*sizeof(double));
  
  u_type->x_Len = vessel_num;
  u_type->y_Len = vessel_num;
  u_type->z_Len = vessel_num;
  u_type->x = malloc(u_type->x_Len*sizeof(double));
  u_type->y = malloc(u_type->y_Len*sizeof(double));
  u_type->z = malloc(u_type->z_Len*sizeof(double));

  y_type->Fx_Len = vessel_num;
  y_type->Fy_Len = vessel_num;
  y_type->Fz_Len = vessel_num;
  y_type->Fx = malloc(y_type->Fx_Len*sizeof(double));
  y_type->Fy = malloc(y_type->Fy_Len*sizeof(double));
  y_type->Fz = malloc(y_type->Fz_Len*sizeof(double));

  return MAP_SAFE;
};


MAP_ERROR_CODE compare_integer_length(const int a, const int b)
{
  if (a!=b) {
    return MAP_FATAL;
  }; 
  return MAP_SAFE;
};


MAP_ERROR_CODE set_node_list(const MAP_ParameterType_t* p_type,  MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, Domain* domain, struct bstrList* node_input_string, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  int i_parsed = 0;
  int next = 0;
  const int num_nodes = node_input_string->qty;
  int fix_num = 0;
  int vessel_num = 0;
  int connect_num = 0;
  Node new_node;
  Node* node_iter = NULL;
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  bstring alias = NULL;
  bstring value_string = NULL;
  const double depth = p_type->depth;
  ReferencePoint u_reference_point;

  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   

  success = allocate_types_for_nodes(u_type, z_type, other_type, y_type, domain, node_input_string, map_msg, ierr);
  success = reset_node(&new_node); /* create an empty node */
   
  for(i=0 ; i<num_nodes ; i++) {         
    list_append(&domain->node, &new_node); /* append node to list */
    node_iter = (Node*)list_get_at(&domain->node, i);
    // success = set_node_vartype(node_iter);
    i_parsed = 0;
    next = 0;
    parsed = bsplits(node_input_string->entry[i], &tokens);

    MAP_BEGIN_ERROR_LOG;  

    while (i_parsed<parsed->qty-1) { /* iterating through all strings */              
      if (parsed->entry[i_parsed]->slen) { /* if the string length is not 0 */
        if (next==0) {            
          next++;
        } else if (next==1) {
          if (biseqcstrcaseless(parsed->entry[i_parsed],"FIX")) {
            node_iter->type = FIX;
            fix_num++;                   /* VarTypePtr              FAST derived  array index */
            success = associate_vartype_ptr(&node_iter->position_ptr.x, other_type->x, fix_num);
            success = associate_vartype_ptr(&node_iter->position_ptr.y, other_type->y, fix_num);
            success = associate_vartype_ptr(&node_iter->position_ptr.z, other_type->z, fix_num);
            success = associate_vartype_ptr(&node_iter->sum_force_ptr.fx, other_type->Fx_anchor, fix_num);
            success = associate_vartype_ptr(&node_iter->sum_force_ptr.fy, other_type->Fy_anchor, fix_num);
            success = associate_vartype_ptr(&node_iter->sum_force_ptr.fz, other_type->Fz_anchor, fix_num);
          } else if (biseqcstrcaseless(parsed->entry[i_parsed],"CONNECT")) {
            node_iter->type = CONNECT;
            connect_num++;
            success = associate_vartype_ptr(&node_iter->position_ptr.x, z_type->x, connect_num);
            success = associate_vartype_ptr(&node_iter->position_ptr.y, z_type->y, connect_num);
            success = associate_vartype_ptr(&node_iter->position_ptr.z, z_type->z, connect_num);
            success = associate_vartype_ptr(&node_iter->sum_force_ptr.fx, other_type->Fx_connect, connect_num);
            success = associate_vartype_ptr(&node_iter->sum_force_ptr.fy, other_type->Fy_connect, connect_num);
            success = associate_vartype_ptr(&node_iter->sum_force_ptr.fz, other_type->Fz_connect, connect_num);
          } else if (biseqcstrcaseless(parsed->entry[i_parsed],"VESSEL")) {
            node_iter->type = VESSEL;
            vessel_num++;
            u_reference_point.x = NULL;
            u_reference_point.y = NULL;
            u_reference_point.z = NULL;
            success = associate_vartype_ptr(&node_iter->position_ptr.x, u_type->x, vessel_num);
            success = associate_vartype_ptr(&node_iter->position_ptr.y, u_type->y, vessel_num);
            success = associate_vartype_ptr(&node_iter->position_ptr.z, u_type->z, vessel_num);                           
            success = associate_vartype_ptr(&node_iter->sum_force_ptr.fx, y_type->Fx, vessel_num);
            success = associate_vartype_ptr(&node_iter->sum_force_ptr.fy, y_type->Fy, vessel_num);
            success = associate_vartype_ptr(&node_iter->sum_force_ptr.fz, y_type->Fz, vessel_num);
            u_reference_point.x = &node_iter->position_ptr.x; /* create reference to input type; this is the convenient update point when u is interpolated in FAST */
            u_reference_point.y = &node_iter->position_ptr.y; /* create reference to input type; this is the convenient update point when u is interpolated in FAST */
            u_reference_point.z = &node_iter->position_ptr.z; /* create reference to input type; this is the convenient update point when u is interpolated in FAST */
            list_append(&domain->u_update_list, &u_reference_point); /* push onto the update list */
          } else {
            set_universal_error_with_message(map_msg, ierr, MAP_FATAL_25, "Value: <%s>", parsed->entry[i_parsed]->data);
          };
          next++;
        } else if (next==2) { /* set initial X node position values */
          alias = bformat("X[%d]", i+1);                          
          success = set_vartype_ptr("[m]", alias, i, &node_iter->position_ptr.x, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_17);
          bdestroy(alias);
          next++;
        } else if (next==3) { /* set initial Y node position values */
          alias = bformat("Y[%d]", i+1);                          
          success = set_vartype_ptr("[m]", alias, i, &node_iter->position_ptr.y, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_18);
          bdestroy(alias);
          next++;
        } else if (next==4) { /* set initial Z node position values */
          alias = bformat("Z[%d]", i+1);                          
          if (biseqcstrcaseless(parsed->entry[i_parsed],"DEPTH")) {         
            if (node_iter->type!=FIX) { /* can only use 'DEPTH' flag in input file for FIX (anchor) nodes */
              set_universal_error_with_message(map_msg, ierr, MAP_FATAL_71, "Value: <%s>", parsed->entry[i_parsed]->data);
            } else {
              value_string = bformat("%f", -depth);                          
              success = set_vartype_ptr("[m]", alias, i, &node_iter->position_ptr.z, value_string); CHECKERRQ(MAP_FATAL_19);
              success = bdestroy(value_string);
            };
          } else { /* all other nodes not using the 'DEPTH' flag */
            success = set_vartype_ptr("[m]", alias, i, &node_iter->position_ptr.z, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_19);
          };        
          bdestroy(alias);
          next++;
        } else if (next==5) { /* set the node mass */            
          alias = bformat("M[%d]", i+1);                          
          success = set_vartype("[kg]", alias, i, &node_iter->M_applied, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_20);
          bdestroy(alias);
          next++;  
        } else if (next==6) { /* set the node buoyancy */
          alias = bformat("B[%d]", i+1);                          
          success = set_vartype("[m^3]", alias, i, &node_iter->B_applied, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_21);
          bdestroy(alias);
          next++; 
        } else if (next==7) { /* set applied X external force (or user guess) of the node */                    
          alias = bformat("FX[%d]", i+1);                          
          success = set_vartype("[N]", alias, i, &node_iter->external_force.fx, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_22);
          bdestroy(alias);
          next++;
        } else if (next==8) { /* set applied Y external force (or user guess) of the node */            
          alias = bformat("FY[%d]", i+1);                          
          success = set_vartype("[N]", alias, i, &node_iter->external_force.fy, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_23);
          bdestroy(alias);
          next++;
        } else if (next==9) { /* set applied Z external force (or user guess) of the node */
          alias = bformat("FZ[%d]", i+1);                          
          success = set_vartype("[N]", alias, i, &node_iter->external_force.fz, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_24);
          bdestroy(alias);
          next++;
        } else {            
          next++;
        };
      };
      i_parsed++;
    };
  
    MAP_END_ERROR_LOG;   

    success = bstrListDestroy(parsed);
    /* @todo: need to make sure next==9; otherwise not enough inputs and an error should
     *        be thrown
     */
  };  

  /* check to make sure the number of allocated array spaces for fortran derived types matches 
   * what was actually set in the node initialization front end.
   */     
  MAP_BEGIN_ERROR_LOG;

  success = compare_integer_length(other_type->Fx_connect_Len, connect_num); CHECKERRQ(MAP_FATAL_49);
  success = compare_integer_length(other_type->Fy_connect_Len, connect_num); CHECKERRQ(MAP_FATAL_49);
  success = compare_integer_length(other_type->Fz_connect_Len, connect_num); CHECKERRQ(MAP_FATAL_49);

  success = compare_integer_length(other_type->Fx_anchor_Len, fix_num); CHECKERRQ(MAP_FATAL_49); // @todo: change error code
  success = compare_integer_length(other_type->Fy_anchor_Len, fix_num); CHECKERRQ(MAP_FATAL_49); // @todo: change error code
  success = compare_integer_length(other_type->Fz_anchor_Len, fix_num); CHECKERRQ(MAP_FATAL_49); // @todo: change error code

  success = compare_integer_length(other_type->x_Len, fix_num); CHECKERRQ(MAP_FATAL_49);
  success = compare_integer_length(other_type->y_Len, fix_num); CHECKERRQ(MAP_FATAL_49);
  success = compare_integer_length(other_type->z_Len, fix_num); CHECKERRQ(MAP_FATAL_49);
              
  success = compare_integer_length(u_type->x_Len, vessel_num); CHECKERRQ(MAP_FATAL_50);
  success = compare_integer_length(u_type->y_Len, vessel_num); CHECKERRQ(MAP_FATAL_50);
  success = compare_integer_length(u_type->z_Len, vessel_num); CHECKERRQ(MAP_FATAL_50);    
              
  success = compare_integer_length(y_type->Fx_Len, vessel_num); CHECKERRQ(MAP_FATAL_51);
  success = compare_integer_length(y_type->Fy_Len, vessel_num); CHECKERRQ(MAP_FATAL_51);
  success = compare_integer_length(y_type->Fz_Len, vessel_num); CHECKERRQ(MAP_FATAL_51);    
              
  success = compare_integer_length(z_type->x_Len, connect_num); CHECKERRQ(MAP_FATAL_52);
  success = compare_integer_length(z_type->y_Len, connect_num); CHECKERRQ(MAP_FATAL_52);
  success = compare_integer_length(z_type->z_Len, connect_num); CHECKERRQ(MAP_FATAL_52);    

  MAP_END_ERROR_LOG;  

  MAP_RETURN_STATUS(*ierr);
};


MAP_ERROR_CODE set_vartype_float(const char* unit, const char* alias, const int num, VarType* type, const double value)
{
  type->name = bfromcstr(alias);
  type->units = bfromcstr(unit); 
  type->ref_counter = 0;
  type->id = num;
  type->value = value;

  return MAP_SAFE;
};


MAP_ERROR_CODE set_vartype(const char* unit, bstring alias, const int num, VarType* type, bstring property)
{
  type->name = bstrcpy(alias);
  type->units = bfromcstr(unit); 
  type->ref_counter = 0;
  type->id = num;
  
  if (!property) { /* this option should only be called for setting line vartypes */
    type->value = -999.9;
  } else {
    if (property->data[0]=='#') { /* this variable is an iterated parameter */      
      type->is_fixed = false;
      if (property->slen==1) { /* implies that property->data = "#" */
        type->value = -999.9;
      } else if (is_numeric(remove_first_character((char*)property->data))) {
        type->value = (double)atof(remove_first_character((char*)property->data));
        type->user_initial_guess = true;
      } else {
        return MAP_FATAL;
      };
    } else { /* this variable is constant */    
      type->is_fixed = true;
      if (is_numeric((char*)property->data)) { 
        type->value = (double)atof((char*)property->data);
      } else {
        return MAP_FATAL;
      };
    };
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE set_vartype_ptr(const char* unit, bstring alias, const int num, VarTypePtr* type, bstring property)
{
  type->name = bstrcpy(alias);
  type->units = bfromcstr(unit); 
  type->ref_counter = 0;
  type->id = num;  

  if (property->data[0]=='#') { /* this variable is an iterated parameter */      
    type->is_fixed = false;
    if (property->slen==1) { /* implies that property->data = "#" */
      *type->value = -999.9;
    } else if (is_numeric(remove_first_character((char*)property->data))) { 
      *type->value = (double)atof(remove_first_character((char*)property->data));
    } else {
      return MAP_FATAL;
    };
  } else { /* this variable is constant */    
    type->is_fixed = true;
    if ((char*)is_numeric((char*)property->data)) { 
      *type->value = (double)atof((char*)property->data);
    } else {
      return MAP_FATAL;
    };
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE set_line_option_flags(struct bstrList* words, int* i_parsed, Line* line_ptr, char* map_msg, MAP_ERROR_CODE* ierr)
{
  // MAP_ERROR_CODE success = MAP_SAFE;
  int index = *i_parsed;
  
  if (biseqcstrcaseless(words->entry[index], "GX_POS")) {
    line_ptr->options.gx_pos_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "GY_POS")) {
    line_ptr->options.gy_pos_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "GZ_POS")) {
    line_ptr->options.gz_pos_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "GX_A_POS")) {
    line_ptr->options.gx_anchor_pos_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "GY_A_POS")) {
    line_ptr->options.gy_anchor_pos_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "GZ_A_POS")) {
    line_ptr->options.gz_anchor_pos_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "GX_FORCE")) {
    line_ptr->options.gx_force_flag= true;
  } else if (biseqcstrcaseless(words->entry[index], "GY_FORCE")) {
    line_ptr->options.gy_force_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "GZ_FORCE")) {
    line_ptr->options.gz_force_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "H_FAIR")) {
    line_ptr->options.H_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "H_ANCH")) {
    line_ptr->options.H_anchor_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "V_FAIR")) {
    line_ptr->options.V_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "V_ANCH")) {
    line_ptr->options.V_anchor_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "TENSION_FAIR")) {
    line_ptr->options.fairlead_tension_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "TENSION_ANCH")) {
    line_ptr->options.anchor_tension_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "X_EXCURSION")) {
    line_ptr->options.horizontal_excursion_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "Z_EXCURSION")) {
    line_ptr->options.vertical_excursion_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "AZIMUTH")) {
    line_ptr->options.azimuth_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "ALTITUDE")) {
    line_ptr->options.altitude_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "ALTITUDE_ANCH")) {
    line_ptr->options.altitude_anchor_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "LINE_TENSION")) {
    line_ptr->options.line_tension_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "OMIT_CONTACT")) {
    line_ptr->options.omit_contact = true;
  } else if (biseqcstrcaseless(words->entry[index], "LINEAR_SPRING")) {
    line_ptr->options.linear_spring = true;
  } else if (biseqcstrcaseless(words->entry[index], "SEG_SIZE")) {
    do {
      index++;
      if (words->qty-1<=index) {
        break;
      };
    } while (words->entry[index]->slen<1);
    if (is_numeric((char*)words->entry[index]->data)) {
      line_ptr->segment_size = (int)atoi((char*)words->entry[index]->data);
      *i_parsed = index;
    } else { /* should not cancel the simulation; simply ignore it */      
      set_universal_error_with_message(map_msg, ierr, MAP_FATAL_18, "Option <%s>", words->entry[index]->data);
    };          
  } else if (biseqcstrcaseless(words->entry[index], "LAY_LENGTH")) {
    line_ptr->options.lay_length_flag = true;
  } else if (biseqcstrcaseless(words->entry[index], "DAMAGE_TIME")) {
    do {
      index++;
      if (words->qty-1<=index) {
        break;
      };
    } while (words->entry[index]->slen<1);
    if (is_numeric((char*)words->entry[index]->data)) {
      line_ptr->options.damage_time_flag = true;
      line_ptr->damage_time = (double)atof((char*)words->entry[index]->data);
      *i_parsed = index;
    } else { /* should not cancel the simulation; simply ignore it */      
      set_universal_error_with_message(map_msg, ierr, MAP_ERROR_1, "Option <%s>", words->entry[index]->data);
    };          
  } else if (biseqcstrcaseless(words->entry[index], "DIAGNOSTIC")) {
    do {
      index++;
      if (words->qty-1<=index) {
        break;
      };
    } while (words->entry[index]->slen<1);
    if (is_numeric((char*)words->entry[index]->data)) {
      line_ptr->options.diagnostics_flag = true;
      line_ptr->diagnostic_type = (int)atoi((char*)words->entry[index]->data);
      *i_parsed = index;
    } else { /* should not cancel the simulation; simply ignore it */      
      set_universal_error_with_message(map_msg, ierr, MAP_ERROR_14, "Option <%s>", words->entry[index]->data);
    };          
  } else {
    /* should not cancel the simulation; simply ignore it */
    set_universal_error_with_message(map_msg, ierr, MAP_WARNING_3, "Option <%s>", words->entry[index]->data);
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE set_line_list(MAP_ConstraintStateType_t* z_type, Domain* domain, struct bstrList* line_input_string, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  int i_parsed = 0;
  int next = 0;
  const int num_lines = line_input_string->qty;
  Line new_line;
  Line* line_iter = NULL;
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  bstring alias = NULL;

  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   
  success = reset_line(&new_line);

  z_type->H_Len = num_lines; 
  z_type->V_Len = num_lines; 
  z_type->H = malloc(z_type->H_Len*sizeof(double));
  z_type->V = malloc(z_type->V_Len*sizeof(double));

  if (z_type->H==NULL || z_type->V==NULL) {
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_53, "Failed allocation of a z_type");
    return MAP_FATAL;
  };
  
  for(i=0 ; i<num_lines ; i++) {         
    list_append(&domain->line, &new_line);
    line_iter = (Line*)list_get_at(&domain->line, i);
    // success = set_line_vartype(line_iter, i); /* @todo: check error */ @rm

    i_parsed = 0;
    next = 0;
    parsed = bsplits(line_input_string->entry[i], &tokens);

    MAP_BEGIN_ERROR_LOG;  

    while (i_parsed<parsed->qty-1) { /* iterating through all strings */              
      if (parsed->entry[i_parsed]->slen) { /* if the string length is not 0 */
        if (next==0) { /* use this first option as an opportunity to set the run-time flags to false */             
          success = associate_vartype_ptr(&line_iter->H, z_type->H, i+1);
          success = associate_vartype_ptr(&line_iter->V, z_type->V, i+1);

          line_iter->H.is_fixed = false;
          alias = bformat("H[%d]", i+1);
          success = set_vartype_ptr("[N]", alias, i, &line_iter->H, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_32);            
          success = bdestroy(alias);

          line_iter->V.is_fixed = false;
          alias = bformat("V[%d]", i+1);
          success = set_vartype_ptr("[N]", alias, i, &line_iter->V, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_32);                        
          success = bdestroy(alias);
             
          next++;
        } else if (next==1) {
          success = associate_line_with_cable_property(line_iter, domain, (char*)parsed->entry[i_parsed]->data, map_msg, ierr); CHECKERRQ(MAP_FATAL_32);           
          next++;
        } else if (next==2) { 
          alias = bformat("Lu[%d]", i+1);
          success = set_vartype("[m]", alias, i, &line_iter->Lu, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_26);
          success = bdestroy(alias);
          next++;
        } else if (next==3) { 
          success = associate_line_with_anchor_node(line_iter, domain, i+1, (char*)parsed->entry[i_parsed]->data,  map_msg, ierr); CHECKERRQ(MAP_FATAL_32);        
          next++;
        } else if (next==4) { 
          success = associate_line_with_fairlead_node(line_iter, domain, i+1, (char*)parsed->entry[i_parsed]->data,  map_msg, ierr); CHECKERRQ(MAP_FATAL_32);        
          next++;
        } else { /* set the node mass */            
          success = set_line_option_flags(parsed, &i_parsed, line_iter, map_msg, ierr);
        };
      };
      i_parsed++;
    };

    MAP_END_ERROR_LOG;   

    success = bstrListDestroy(parsed);
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE push_variable_to_output_list(OutputList* y_list, const int i, double* variable_ref, const char* alias, const char* units) {
  int size = 0;

  VarTypePtr vartype_ptr;
  VarTypePtr* iter_vartype = NULL;

  list_append(&y_list->out_list_ptr, &vartype_ptr); /* @todo: this is not correct. Should point to fairlead->sumForce.fy */
  size = list_size(&y_list->out_list_ptr);
  iter_vartype = (VarTypePtr*)list_get_at(&y_list->out_list_ptr, size-1);
  iter_vartype->value = variable_ref;
  iter_vartype->name = bformat("%s[%d]", alias, i);
  iter_vartype->units = bformat("%s", units);      

  return MAP_SAFE;
};


MAP_ERROR_CODE set_output_list(Domain* domain, MAP_InitOutputType_t* io_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE  success = MAP_SAFE;
  Line* line_iter = NULL;
  OutputList* y_list = domain->y_list;
  // int size = 0;
  int line_num = 1;
  // VarTypePtr* iter_vartype = NULL;

  list_iterator_start(&domain->line); /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->line)) { /* tell whether more values available */ 
    line_iter = (Line*)list_iterator_next(&domain->line);    
    
    if (line_iter->options.gx_anchor_pos_flag) {
      list_append(&y_list->out_list_ptr, &line_iter->anchor->position_ptr.x);      
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.gy_anchor_pos_flag) {
      list_append(&y_list->out_list_ptr, &line_iter->anchor->position_ptr.y);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.gz_anchor_pos_flag) {
      list_append(&y_list->out_list_ptr, &line_iter->anchor->position_ptr.z);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.gx_pos_flag) {
      list_append(&y_list->out_list_ptr, &line_iter->fairlead->position_ptr.x);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.gy_pos_flag) {
      list_append(&y_list->out_list_ptr, &line_iter->fairlead->position_ptr.y);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.gz_pos_flag) {
      list_append(&y_list->out_list_ptr, &line_iter->fairlead->position_ptr.z);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.H_flag) {
      list_append(&y_list->out_list_ptr, &line_iter->H);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.V_flag) {
      list_append(&y_list->out_list_ptr, &line_iter->V);
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };
    
    if (line_iter->options.H_anchor_flag) {
      success = push_variable_to_output_list(y_list, line_num, &line_iter->H_at_anchor, "H_a", "[N]");
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.V_anchor_flag) {
      success = push_variable_to_output_list(y_list, line_num, &line_iter->V_at_anchor, "V_a", "[N]");
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.gx_force_flag) {      
      success = push_variable_to_output_list(y_list, line_num, &line_iter->fx_fairlead, "Fx", "[N]");
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.gy_force_flag) {
      success = push_variable_to_output_list(y_list, line_num, &line_iter->fy_fairlead, "Fy", "[N]");
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.gz_force_flag) {
      success = push_variable_to_output_list(y_list, line_num, &line_iter->fz_fairlead, "Fz", "[N]");
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.fairlead_tension_flag) {
      success = push_variable_to_output_list(y_list, line_num, &line_iter->T, "T", "[N]");
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.anchor_tension_flag) {
      success = push_variable_to_output_list(y_list, line_num, &line_iter->T_at_anchor, "T_a", "[N]");
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.lay_length_flag) {
      success = push_variable_to_output_list(y_list, line_num, &line_iter->Lb, "Lb", "[m]");
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.horizontal_excursion_flag) {
      success = push_variable_to_output_list(y_list, line_num, &line_iter->l, "l", "[m]");
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.vertical_excursion_flag) {
      success = push_variable_to_output_list(y_list, line_num, &line_iter->h, "h", "[m]");
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.azimuth_flag) {
      success = push_variable_to_output_list(y_list, line_num, &line_iter->psi, "psi", "[rad]");
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.altitude_flag) {
      success = push_variable_to_output_list(y_list, line_num, &line_iter->alpha, "alpha", "[rad]");
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    if (line_iter->options.altitude_anchor_flag) {
      success = push_variable_to_output_list(y_list, line_num, &line_iter->alpha_at_anchor, "alpha_a", "[rad]");
      io_type->writeOutputHdr_Len++;
      io_type->writeOutputUnt_Len++;
    };

    line_num++;
  };
  list_iterator_stop(&domain->line); /* ending the iteration session */  

  return MAP_SAFE;
};


MAP_ERROR_CODE reset_line(Line* line_ptr)
{
  /* run-time flags */
  line_ptr->options.gx_pos_flag = false;
  line_ptr->options.gy_pos_flag = false;
  line_ptr->options.gz_pos_flag = false;
  line_ptr->options.gx_anchor_pos_flag = false;
  line_ptr->options.gy_anchor_pos_flag = false;
  line_ptr->options.gz_anchor_pos_flag = false;
  line_ptr->options.gx_force_flag = false;
  line_ptr->options.gy_force_flag = false;
  line_ptr->options.gz_force_flag = false;
  line_ptr->options.H_flag = false;
  line_ptr->options.H_anchor_flag = false;
  line_ptr->options.V_flag = false;
  line_ptr->options.V_anchor_flag = false;
  line_ptr->options.fairlead_tension_flag = false;
  line_ptr->options.anchor_tension_flag = false;
  line_ptr->options.horizontal_excursion_flag = false;
  line_ptr->options.vertical_excursion_flag = false;
  line_ptr->options.azimuth_flag = false;
  line_ptr->options.altitude_flag = false;
  line_ptr->options.altitude_anchor_flag = false;
  line_ptr->options.line_tension_flag = false;
  line_ptr->options.omit_contact = false;
  line_ptr->options.lay_length_flag = false;
  line_ptr->options.damage_time_flag = false;
  line_ptr->options.diagnostics_flag = false;
  line_ptr->options.linear_spring = false;
  
  line_ptr->line_property = NULL;      
  line_ptr->label = NULL;
  line_ptr->line_tension  = NULL;
  line_ptr->anchor = NULL;             /* Anchor node */
  line_ptr->fairlead = NULL;           /* Fairlead node */
  
  line_ptr->H.name = NULL;
  line_ptr->H.units = NULL;
  line_ptr->V.name = NULL;
  line_ptr->V.units = NULL;

  line_ptr->psi = -999.9;
  line_ptr->alpha = -999.9;
  line_ptr->alpha_at_anchor = -999.9;
  line_ptr->l = -999.9;
  line_ptr->Lb = -999.9;
  line_ptr->h = -999.9;
  line_ptr->H.value = NULL;
  line_ptr->V.value = NULL;
  line_ptr->H_at_anchor = -999.9;
  line_ptr->V_at_anchor = -999.9;
  line_ptr->fx_fairlead = -999.9;
  line_ptr->fy_fairlead = -999.9;
  line_ptr->fz_fairlead = -999.9;
  line_ptr->fx_anchor = -999.9;
  line_ptr->fy_anchor = -999.9;
  line_ptr->fz_anchor = -999.9;
  line_ptr->T = -999.9;
  line_ptr->T_at_anchor = -999.9;

  line_ptr->residual_norm = 999.9;
  line_ptr->damage_time = -999.9;
  line_ptr->diagnostic_type = -9999;
  line_ptr->segment_size = 10;
  return MAP_SAFE;
};


MAP_ERROR_CODE reset_node(Node* node_ptr)
{
  node_ptr->position_ptr.x.name = NULL;
  node_ptr->position_ptr.x.units = NULL;
  node_ptr->position_ptr.x.value = NULL;
  node_ptr->position_ptr.y.name = NULL;
  node_ptr->position_ptr.y.units = NULL;
  node_ptr->position_ptr.y.value = NULL;
  node_ptr->position_ptr.z.name = NULL;
  node_ptr->position_ptr.z.units = NULL;
  node_ptr->position_ptr.z.value = NULL;
  node_ptr->M_applied.name = NULL;
  node_ptr->M_applied.units = NULL;
  node_ptr->B_applied.name = NULL;
  node_ptr->B_applied.units = NULL;
  node_ptr->sum_force_ptr.fx.name = NULL;
  node_ptr->sum_force_ptr.fx.units = NULL;
  node_ptr->sum_force_ptr.fx.value = NULL;
  node_ptr->sum_force_ptr.fy.name = NULL;
  node_ptr->sum_force_ptr.fy.units = NULL;
  node_ptr->sum_force_ptr.fy.value = NULL;
  node_ptr->sum_force_ptr.fz.name = NULL;
  node_ptr->sum_force_ptr.fz.units = NULL;
  node_ptr->sum_force_ptr.fz.value = NULL;
 
  node_ptr->external_force.fx.name = NULL;
  node_ptr->external_force.fx.units = NULL;
  node_ptr->external_force.fy.name = NULL;
  node_ptr->external_force.fy.units = NULL;
  node_ptr->external_force.fz.name = NULL;
  node_ptr->external_force.fz.units = NULL;

  node_ptr->M_applied.value = -999.9;
  node_ptr->B_applied.value = -999.9; 
  node_ptr->external_force.fx.value = -999.9;
  node_ptr->external_force.fy.value = -999.9;
  node_ptr->external_force.fz.value = -999.9;
  node_ptr->external_force.fx.user_initial_guess = false;
  node_ptr->external_force.fy.user_initial_guess = false;
  node_ptr->external_force.fz.user_initial_guess = false;

  node_ptr->sum_force_ptr.fx.is_fixed = false;
  node_ptr->sum_force_ptr.fy.is_fixed = false;
  node_ptr->sum_force_ptr.fz.is_fixed = false;
  return MAP_SAFE;
};


MAP_ERROR_CODE associate_line_with_cable_property(Line* line_ptr, Domain* domain, const char* word, char* map_msg, MAP_ERROR_CODE* ierr)
{
  // MAP_ERROR_CODE success = MAP_SAFE;
  CableLibrary* library_iterator = NULL;

  library_iterator = NULL;
  line_ptr->line_property = NULL;

  list_iterator_start(&domain->library); /* starting an iteration session */
  while (list_iterator_hasnext(&domain->library)) { /* tell whether more values available */
    library_iterator = (CableLibrary*)list_iterator_next(&domain->library);
    if (biseqcstrcaseless(library_iterator->label, word)) {      
      line_ptr->line_property = library_iterator;
      list_iterator_stop(&domain->library); /* ending the iteration session */  
      break;
    }; 
  };
  list_iterator_stop(&domain->library); /* ending the iteration session */  
  if (line_ptr->line_property==NULL) {        
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_27, "No libraries match <%s>.", word);
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE associate_line_with_anchor_node(Line* line_ptr, Domain* domain, const int line_num, const char* word, char* map_msg, MAP_ERROR_CODE* ierr)
{
  // MAP_ERROR_CODE success = MAP_SAFE;
  Node* node_iter = NULL;
  int node_num = 0;

  line_ptr->anchor = NULL;
  
  if (is_numeric(word)) {
    node_num = (int)atoi(word); 
    node_iter = (Node*)list_get_at(&domain->node, node_num-1);
    line_ptr->anchor = node_iter; /* create the associate with anchor here */
    if (!node_iter) {
      set_universal_error_with_message(map_msg, ierr, MAP_FATAL_30, "Line %d.", line_num);
      return MAP_FATAL;
    };
  } else {
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_28, "Line %d.", line_num);
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE associate_line_with_fairlead_node(Line* line_ptr, Domain* domain, const int line_num, const char* word, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Node* node_iter = NULL;
  int node_num = 0;
  // MAP_ERROR_CODE success = MAP_SAFE;

  line_ptr->fairlead = NULL;

  if (is_numeric(word)) {
    node_num = (int)atoi(word); 
    node_iter = (Node*)list_get_at(&domain->node, node_num-1);
    line_ptr->fairlead = node_iter; /* create the associate with anchor here */
    if (!node_iter) {
      set_universal_error_with_message(map_msg, ierr, MAP_FATAL_31, "Line %d.", line_num);
      return MAP_FATAL;
    };
  } else {
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_29, "Line %d.", line_num);
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE is_numeric(const char* string)
{
  char* p = NULL;
  if (string==NULL || *string=='\0' || isspace(*string)) {
    return MAP_SAFE;
  };
  strtod (string, &p);
  if (*p=='\0') {
    return MAP_FATAL;
  } else {
    return MAP_SAFE;
  };
};


void log_initialization_information(MAP_InitInputType_t* init_type, MAP_ParameterType_t* p_type, MAP_OutputType_t* y_type, MAP_OtherStateType_t* other_type, Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  InitializationData* init_data = init_type->object;   

  MAP_BEGIN_ERROR_LOG; 
  if ( init_data->summary_file_name->data[0] ) { // don't write this file if the file name isn't specified
    success = write_summary_file(init_data, p_type, domain, map_msg, ierr); CHECKERRQ(MAP_FATAL_37); 
  }

  success = write_summary_file(init_data, p_type, domain, map_msg, ierr); CHECKERRQ(MAP_FATAL_37);           
  success = get_iteration_output_stream(y_type, other_type, map_msg, ierr); // @todo CHECKERRQ()    
  MAP_END_ERROR_LOG; 
};


MAP_ERROR_CODE associate_vartype_ptr(VarTypePtr* type, double* arr, int index)
{
  type->value = &arr[index-1];
  return MAP_SAFE;
};


void copy_target_string(char* target, unsigned char* source)
{
  while (*source) {
    *target = *source;
    source++;
    target++;
  };
  *target = '\0';
};


MAP_ERROR_CODE map_get_version(MAP_InitOutputType_t* io_type)
{
  bstring out_string = NULL;
  int ret = 0;

  /* first set the program version defined in the mapsys.h header file 
   * @todo: program version should be tied to the gi revision numner
   */
  out_string = bformat("<%s>",PROGVERSION);
  if (out_string->slen>MAX_INIT_VERSION_STRING_LENGTH) { /* overflow */
    return MAP_FATAL; /* @todo: give proper error code */
  };
  copy_target_string(io_type->version, out_string->data);
  ret = bdestroy(out_string);

  /* the set the compiling date. This is #defined in the mapsys.h header */
  out_string = bformat("<%c%c%c-%c%c-%c%c%c%c>",BUILD_MONTH_CH0,BUILD_MONTH_CH1,BUILD_MONTH_CH2,BUILD_DAY_CH0,BUILD_DAY_CH1,BUILD_YEAR_CH0,BUILD_YEAR_CH1,BUILD_YEAR_CH2,BUILD_YEAR_CH3);
  if (out_string->slen>MAX_INIT_COMPILING_DATA_STRING_LENGTH) { /* overflow */
    return MAP_FATAL; /* @todo: give proper error code */
  };
  copy_target_string(io_type->compilingData, out_string->data);
  ret = bdestroy(out_string);
  return MAP_SAFE;
};


void print_machine_name_to_screen( ) {
  // __get_machine_name(name);
  printf( "%s Ver. %s ", PROGNAME, PROGVERSION); 
  printf( "%c",BUILD_MONTH_CH0 );// build month
  printf( "%c",BUILD_MONTH_CH1 );
  printf( "%c",BUILD_MONTH_CH2 );
  printf( "-" );
  printf( "%c",BUILD_DAY_CH0 );// build day
  printf( "%c",BUILD_DAY_CH1 );
  printf( "-" );
  printf( "%c",BUILD_YEAR_CH0 ); // build year 
  printf( "%c",BUILD_YEAR_CH1 );
  printf( "%c",BUILD_YEAR_CH2 );
  printf( "%c\n",BUILD_YEAR_CH3 );
}


const char* remove_first_character(const char* string)
{
  return string+1;
};


MAP_ERROR_CODE print_help_to_screen()
{
  print_machine_name_to_screen( );

  printf("MAP Input file section definitions:\n");
  printf("  Line dictionary definitions:\n");   
  printf("    -LineType, --User-defined name of line [-]  \n");   
  printf("    -Diam,     --Line diameter, used to calculate area and line displacement per unit length [m]  \n");   
  printf("    -MassDen,  --Mass (in air) per unit length [kg/m]  \n");   
  printf("    -EA,       --Axial stiffness [N] \n");   
  printf("    -CB,       --Cable/seabed Coulumb friction coefficient [-]  \n");   
  printf("    -CIntDamp, --Internal structural damping coefficient [Pa-s]  \n");   
  printf("    -Ca,       --Cross-flow added-mass coefficient [-]\n");   
  printf("    -Cdn,      --Cross-flow drag coefficient [-]\n");   
  printf("    -Cdt,      --Tangent (skin) drag coefficient[-]\n");   
  printf("  Node property definitions:\n");
  printf("    -Node,     --Node number; first starts at 1 [-]\n");
  printf("    -Type,     --Type of node. Must be one of: VESSEL, FIX, CONNECT [-]\n");
  printf("    -X,        --Node X position. '#' must prefix CONNECT nodes; constitutes user initial guess [m]\n");
  printf("    -Y,        --Node Y position. '#' must prefix CONNECT nodes; constitutes user initial guess [m]\n");
  printf("    -Z,        --Node Z position. '#' must prefix CONNECT nodes; constitutes user initial guess [m]\n");
  printf("    -M,        --Applied point mass at node [kg]\n");
  printf("    -B,        --Applied point buoyancy module at node [m^3]\n");  
  printf("    -FX,       --Applied X external force at node. '#' must prefix VESSEL and FIX nodes [N]\n");
  printf("    -FY,       --Applied Y external force at node. '#' must prefix VESSEL and FIX nodes [N]\n");
  printf("    -FZ,       --Applied Z external force at node. '#' must prefix VESSEL and FIX nodes [N]\n");
  printf("  Line property definitions:\n");
  printf("    -Line,     --Line number; first starts at 1 [-]\n");
  printf("    -LineType, --Must match property defined in 'Line Dictions'[-]\n");
  printf("    -UnstrLen, --Unstretched line length [m]\n");
  printf("    -NodeAnch, --Anchor node number corresponding to 'Node Property Definitions' section [-]\n");
  printf("    -NodeFair, --Fairlead node number corresponding to 'Node Property Definitions' section [-]\n");
  printf("    -Flags,    --User run-time flag; see below [-]\n");
  printf("    \n");
  printf("  Line run-time options definitions\n");
  printf("    Outputs:\n");
  printf("      -gx_pos,       --Fairlead position in global X [m]\n");
  printf("      -gy_pos,       --Fairlead position in global Y [m]\n");
  printf("      -gx_pos,       --Fairlead position in global Z [m]\n");
  printf("      -gx_a_pos,     --Anchor position in global X [m]\n");
  printf("      -gy_a_pos,     --Anchor position in global Y [m]\n");
  printf("      -gz_a_pos,     --Anchor position in global Z [m]\n");
  printf("      -gx_force,     --Fairlead force in global X (include applied forces) [N]\n");
  printf("      -gy_force,     --Fairlead force in global Y (include applied forces) [N]\n");
  printf("      -gz_force,     --Fairlead force in global Z (include applied forces) [N]\n");
  printf("      -h_fair,       --Horizontal force at fairlead (does NOT include applied forces) [N]\n");
  printf("      -v_fair,       --Vertical force at fairlead (does NOT include applied forces) [N]\n");
  printf("      -h_anch,       --Horizontal force at anchor (does NOT include applied forces) [N]\n");
  printf("      -v_anch,       --Vertical force at anchor (does NOT include applied forces) [N]\n");
  printf("      -tension_fair, --Line force-magnitude at fairlead (include applied loads) [N]\n");
  printf("      -tension_anch, --Line force-magnitude at anchor (include applied loads) [N]\n");
  printf("      -azimuth,      --Line lateral offset angle global X axis [rad]\n");
  printf("      -altitude,     --Line inclination angle relative to global XY plane at fairlead [rad]\n");
  printf("      -lay_length,   --Length of line on seabed [m]\n");
  printf("      -line_tension, -- \n");
  printf("    Model features:\n");
  printf("      -omit_contact,       --Ignore cable/seabed contact\n");
  printf("      -seg_size <10>,      --Number of discrete lines in line\n");
  printf("      -damage_time <NULL>, --Line breakage occurs at specified time [s]\n");
  printf("      -diagnostic <TIME>,  --Run line solver diagnostics until specified time [s] is reached\n");
  printf("\n");
  printf("  Model option definitions\n");
  printf("    General model features:\n");
  printf("      -ref_position <0.0> <0.0> <0.0>\n");
  printf("      -repeat <NULL> ... <NULL>\n");
  printf("    MSQS solver options:\n");
  printf("      -inner_ftol <float>,\n");
  printf("      -inner_gtol <float>,\n");
  printf("      -inner_xtol <float>,\n");
  printf("      -inner_max_its <int>,\n");
  printf("      -outer_tol <float>,\n");
  printf("      -outer_max_its <int>,\n");
  printf("      -outer_epsilon <float>,\n");
  printf("      -outer_bd,\n");
  printf("      -outer_cd,\n");
  printf("      -outer_fd,\n");
  printf("      -pg_cooked <1000.0> <1.0>,\n");
  printf("      -krylov_accelerator <3>,\n");
  printf("      -integration_dt <0.01>,\n");
  printf("    LM model feature (not suported yet):\n");
  printf("      -kb_default      --Seabed stiffness parameter\n");
  printf("      -cb_default      --Seabed damping parameter\n");
  printf("      -wave_kinematics --Enables wave kinematics to drag interaction from surface waves\n");
  printf("      -lm_model        --Enable the lumped-mass model\n");
  printf( "\nMAP++ Copyright (C) 2014 and Apache'd by Marco Masciola and others\n" );
  printf( "SimCList Copyright (C) 2010 by Mij <http://mij.oltrelinux.com/devel/simclist/>\n" );
  printf( "MinPack Copyright (C) 1999 by the University of Chicago\n" );
  printf( "Modifications to MinPack by Frederic Devernay <http://devernay.free.fr/hacks/cminpack/>\n" );
  printf( "\nMAP++ is free software; see the source for copying conditions.\n" );
  printf( "This software is distributed on an \"AS IS\" BASIS, WITHOUT WARRANTIES\n" );
  printf( "OR CONDITIONS OF ANY KIND, either express or implied. See\n" );
  printf( "<http://www.apache.org/licenses/LICENSE-2.0> forr more details.\n" );
  printf("    \n");
  return MAP_SAFE;
};


MAP_ERROR_CODE associate_constraint_states(Domain* domain, MAP_ConstraintStateType_t* z_type)
{
  Node* node_iter = NULL;
  Line* line_iter = NULL;
  int next = 0;
  MAP_ERROR_CODE success = MAP_SAFE;
  
  MAP_BEGIN_ERROR_LOG;

  list_iterator_start(&domain->line);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->line)) { /* tell whether more values available */
    line_iter = (Line*)list_iterator_next(&domain->line);
    success = associate_vartype_ptr(&line_iter->H, z_type->H, next+1);
    success = associate_vartype_ptr(&line_iter->V, z_type->V, next+1);
    next++;
  };
  list_iterator_stop(&domain->line); /* ending the iteration "session" */    
  
  next = 0;
  list_iterator_start(&domain->node);            /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->node)) { /* tell whether more values available */
    node_iter = (Node*)list_iterator_next(&domain->node);
    if (node_iter->type==CONNECT) {
      success = associate_vartype_ptr(&node_iter->position_ptr.x, z_type->x, next+1);
      success = associate_vartype_ptr(&node_iter->position_ptr.y, z_type->y, next+1);
      success = associate_vartype_ptr(&node_iter->position_ptr.z, z_type->z, next+1);
      next++;
    };
  };
  list_iterator_stop(&domain->node); /* ending the iteration "session" */    

  MAP_END_ERROR_LOG;

  return MAP_SAFE;
}
