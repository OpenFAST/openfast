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


#include "numeric.h"
#include "lineroutines.h"
#include "jacobian.h"
#include "residual.h"


MAP_ERROR_CODE root_finding_step(OuterSolveAttributes* ns, const int n, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, double* error, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  const int z_size = z_type->z_Len; 
  const int THREE = 3;
  int i = 0;

  MAP_BEGIN_ERROR_LOG;
  
  success = lu(ns, n, map_msg, ierr); CHECKERRQ(MAP_FATAL_74);
  success = lu_back_substitution(ns, n, map_msg, ierr); CHECKERRQ(MAP_FATAL_74);
  
  /* Note that: ns->x = J^(-1) * F
   *  [x,y,z]_i+1 =  [x,y,z]_i - J^(-1) * F        
   */   
  for (i=0 ; i<z_size ; i++) { 
    z_type->x[i] -= ns->x[THREE*i];
    z_type->y[i] -= ns->x[THREE*i+1];
    z_type->z[i] -= ns->x[THREE*i+2];
    *error += (pow(other_type->Fx_connect[i],2)+ pow(other_type->Fy_connect[i],2) + pow(other_type->Fz_connect[i],2));
  };

  MAP_END_ERROR_LOG;

  return MAP_SAFE;
};


int inner_function_evals(void* line_ptr, int m, int n, const __cminpack_real__* x, __cminpack_real__* fvec, __cminpack_real__* fjac, int ldfjac, int iflag) 
{
  Line* line = (Line*)line_ptr;
  const double Fh = x[0];
  // const double Fh = fabs(x[0]) > MAP_HORIZONTAL_TOL ? fabs(x[0]) : MAP_HORIZONTAL_TOL;
  const double Fv = x[1];  
  const double EA = line->line_property->EA;
  const double Lu = line->Lu.value;
  const double height = line->h;
  const double length = line->l;
  // const double height = line->h > MAP_HORIZONTAL_TOL ? line->h : MAP_HORIZONTAL_TOL;
  // const double length = line->l > MAP_HORIZONTAL_TOL ? line->l : MAP_HORIZONTAL_TOL;
  const double omega = line->line_property->omega;
  const double cb = line->line_property->cb;
  const bool contactFlag = line->options.omit_contact;
  
  if (iflag==0) {
    return 0;
  };
 
  /* Taken from the preceeding FAST 7 HydroDyn.f90 source (verbatim):
   * To avoid an ill - conditioned situation, ensure that the initial guess for HF is not less than or equal to zero.Similarly, avoid the problems
   * associated with having exactly vertical(so that HF is zero) or exactly horizontal(so that VF is zero) lines by setting the minimum values
   * equal to the tolerance.This prevents us from needing to implement the known limiting solutions for vertical or horizontal lines(and thus
   * complicating this routine):
   * 
   * HF = MAX(HF, Tol)
   * XF = MAX(XF, Tol)
   * ZF = MAX(ZF, TOl)
   */  
  if (iflag!=2) {
    if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      fvec[0] = residual_function_length_no_contact(Fv, Fh, omega, Lu, EA, length);
      fvec[1] = residual_function_height_no_contact(Fv, Fh, omega, Lu, EA, height); 
    } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
      fvec[0] = residual_function_length_contact(Fv, Fh, omega, Lu, EA, length, cb);
      fvec[1] = residual_function_height_contact(Fv, Fh, omega, Lu, EA, height, cb); 
    };
  } else {
    if (contactFlag==true || omega<0.0 || (Fv-omega*Lu)>0.0) { /* true when no portion of the line rests on the seabed */
      fjac[0] = jacobian_dxdh_no_contact(Fv, Fh, omega, Lu, EA);
      fjac[1] = jacobian_dxdv_no_contact(Fv, Fh, omega, Lu, EA);
      fjac[2] = jacobian_dzdh_no_contact(Fv, Fh, omega, Lu, EA);
      fjac[3] = jacobian_dzdv_no_contact(Fv, Fh, omega, Lu, EA);
    } else { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
      fjac[0] = jacobian_dxdh_contact(Fv, Fh, omega, Lu, EA, cb);
      fjac[1] = jacobian_dxdv_contact(Fv, Fh, omega, Lu, EA, cb);
      fjac[2] = jacobian_dzdh_contact(Fv, Fh, omega, Lu, EA, cb);
      fjac[3] = jacobian_dzdv_contact(Fv, Fh, omega, Lu, EA, cb);
    };
  };
  return 0;
};


MAP_ERROR_CODE lu(OuterSolveAttributes* ns, const int n, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int j = 0;
  int k = 0;
           
  for (i=0 ; i<n ; i++) {
    ns->l[i][i] = 1.0;
    for (j=i+1 ; j<n ; j++) {
      if (fabs(ns->jac[i][i])<MACHINE_EPSILON) {
       return MAP_FATAL;
      };
      ns->l[j][i] = (ns->jac[j][i])/(ns->jac[i][i]);
      /* ns->jac[j][j] = ns->l[j][j] */
      for (k=i+1 ; k<n ; k++) {
        ns->jac[j][k] = ns->jac[j][k] - (ns->l[j][i])*(ns->jac[i][k]);
      };
    };
    
    for (k=i ; k<n ; k++) {
      ns->u[i][k] = ns->jac[i][k] ;
    };
  };

  return MAP_SAFE;
};


/**
 * Ax = b -> LUx = b. Then y is defined to be Ux
 */
MAP_ERROR_CODE lu_back_substitution(OuterSolveAttributes* ns, const int n, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int i = 0;
  int j = 0;
  
  /* Forward solve Ly = b */
  for (i=0 ; i<n ; i++) {
    /* @todo: I think y can be eliminated. It is simply a copy of b */
    ns->y[i] = ns->b[i];
    for (j=0 ; j<i ; j++) {
      ns->y[i] -= (ns->l[i][j])*(ns->y[j]);
    };
    if (fabs(ns->l[i][i])<MACHINE_EPSILON) {
      return MAP_FATAL;
    };
    ns->y[i] /= ns->l[i][i];    
  };

  /* Backward solve Ux = y */
  for (i=n-1 ; i>=0 ; i--) {
    ns->x[i] = ns->y[i];
    for (j=i+1 ; j<n ; j++) {
      ns->x[i] -= (ns->u[i][j])*(ns->x[j]);
    };    
    if (fabs(ns->u[i][i])<MACHINE_EPSILON) {
      return MAP_FATAL;
    };
    ns->x[i] /= ns->u[i][i];
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE call_minpack_lmder(Line* line, InnerSolveAttributes* inner_opt, const int line_num, const float time, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;

  /* initial guess vector is set in set_line_initial_guess(..); otherwise, the previous solution is used as the initial guess */
  inner_opt->x[0] = fabs(*(line->H.value)) > MAP_HORIZONTAL_TOL ? fabs(*(line->H.value)) : MAP_HORIZONTAL_TOL;
  // inner_opt->x[0] = *(line->H.value);
  inner_opt->x[1] = *(line->V.value);

  line->evals = 0;
  line->njac_evals = 0;

  inner_opt->info = __cminpack_func__(lmder)(inner_function_evals, 
                                      line, 
                                      inner_opt->m, 
                                      inner_opt->n, 
                                      inner_opt->x, 
                                      inner_opt->fvec, 
                                      inner_opt->fjac, 
                                      inner_opt->ldfjac, 
                                      inner_opt->f_tol, 
                                      inner_opt->x_tol, 
                                      inner_opt->g_tol, 
                                      inner_opt->max_its, 
                                      inner_opt->diag,
                                      inner_opt->mode, 
                                      inner_opt->factor, 
                                      inner_opt->nprint, 
                                      &line->evals, 
                                      &line->njac_evals, 
                                      inner_opt->ipvt, 
                                      inner_opt->qtf, 
                                      inner_opt->wa1 ,
                                      inner_opt->wa2 ,
                                      inner_opt->wa3 , 
                                      inner_opt->wa4);
  
  line->residual_norm = (double)__minpack_func__(enorm)(&inner_opt->m, inner_opt->fvec);
  
  if (line->options.diagnostics_flag && (double)line->diagnostic_type>time /* || line->residual_norm>inner_opt->f_tol */ ) {
    printf("\n      %4.3f [sec]  Line %d\n",time, line_num);
    printf("      ----------------------------------------------------\n");
    printf("      Residual l2 norm at solution:  %15.7g\n", line->residual_norm);
    printf("      Function evaluations:         %10i\n", line->evals);
    printf("      Jacobian evaluations:         %10i\n", line->njac_evals);
	if (line->residual_norm>inner_opt->f_tol) {
		printf("      WARNING: l2 norm is much larger than f_tol. Premature convergence is likely\n");
	}
    printf("      Exit parameter                %10i\n\n", inner_opt->info);
  };
  
  *(line->H.value) = inner_opt->x[0];
  *(line->V.value) = inner_opt->x[1];
  line->converge_reason = inner_opt->info;
  
  switch (inner_opt->info) {
  case 0 :
    success = MAP_FATAL;
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_39, "Line segment %d.", line_num);
    break;
  case 1 :
    success = MAP_SAFE;
    break;
  case 2 :
    success = MAP_SAFE;
    break;
  case 3 :
    success = MAP_SAFE;
    break;
  case 4 :
    success = MAP_SAFE;
    break;
  case 5 :
    success = MAP_FATAL;
    set_universal_error_with_message(map_msg, ierr, MAP_FATAL_40, "Line segment %d.", line_num);
    break;
  case 6 :
    success = MAP_FATAL;
    set_universal_error_with_message(map_msg, ierr, MAP_ERROR_11, "Line segment %d.", line_num);
    break;
  case 7 :
    success = MAP_FATAL;
    set_universal_error_with_message(map_msg, ierr, MAP_ERROR_13, "Line segment %d.", line_num);
    break;
  case 8 :
    success = MAP_FATAL;
    set_universal_error_with_message(map_msg, ierr, MAP_ERROR_12, "Line segment %d.", line_num);
    break;
  default :
    success = MAP_SAFE;
    break;
  };
  return MAP_SAFE;
};
