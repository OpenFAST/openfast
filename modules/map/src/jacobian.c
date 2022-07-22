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


#include "jacobian.h"
#include "lineroutines.h"


double jacobian_dxdh_no_contact(const double V, const double H, const double w, const double Lu, const double EA)
{
  return (ARCSINH(V/H) - ARCSINH((V-w*Lu)/H))/w - ((V/H + pow(V/H, 2)/sqrt(1.0 + pow(V/H, 2)))/(V/H + sqrt(1.0 + pow(V/H, 2))) 
                                                   - ((V-w*Lu)/H + pow((V-w*Lu)/H, 2)/sqrt(1.0 + pow((V-w*Lu)/H, 2)))
                                                   /((V-w*Lu)/H + sqrt(1.0 + pow((V-w*Lu)/H, 2))))/w + (Lu/(EA));      
};


double jacobian_dxdv_no_contact(const double V, const double H, const double w, const double Lu, const double EA)
{
  return ((1.0 + V/H /sqrt(1.0 + pow(V/H, 2)))/(V/H + sqrt(1.0 + pow(V/H, 2))) 
          - (1.0 + (V-w*Lu)/H /sqrt(1.0 + pow( (V-w*Lu)/H , 2)))
          /((V-w*Lu)/H + sqrt(1.0 + pow((V-w*Lu)/H, 2))))/w;
};


double jacobian_dzdh_no_contact(const double V, const double H, const double w, const double Lu, const double EA)
{
  return ( sqrt( 1.0 + pow( V/H , 2) ) - sqrt( 1.0 + pow( (V-w*Lu)/H , 2) ) )/w  
    - ( pow( V/H , 2 )/sqrt( 1.0 + pow( V/H , 2) ) - pow( (V-w*Lu)/H , 2)/sqrt( 1.0 + pow( (V-w*Lu)/H , 2) ) )/w;    
};


double jacobian_dzdv_no_contact(const double V, const double H, const double w, const double Lu, const double EA)
{
  return ( V/H/sqrt( 1.0 + pow( V/H , 2) ) - (V-w*Lu)/H /sqrt( 1.0 + pow( (V-w*Lu)/H , 2) ) )/w + (Lu/(EA));
};


double jacobian_dxdh_contact(const double V, const double H, const double w, const double Lu, const double EA, const double cb)
{
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return log((V/H) + sqrt(1.0 + pow(V/H,2)))/w - (((V/H) + (V/H)*(V/H)/sqrt(1.0 + pow(V/H,2)))/((V/H) + sqrt(1.0 + pow(V/H,2))))/w + (Lu/EA);
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return log((V/H) + sqrt(1.0 + pow(V/H,2)))/w - (((V/H) + (V/H)*(V/H)/sqrt(1.0 + pow(V/H,2)))/((V/H) + sqrt(1.0 + pow(V/H,2))))/w + (Lu/EA) - ((Lu-V/w) - (H/w)/cb)/EA;
  };
};


double jacobian_dxdv_contact(const double V, const double H, const double w, const double Lu, const double EA, const double cb)
{
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return ((1.0 + (V/H)/sqrt(1.0 + pow(V/H,2)))/((V/H) + sqrt(1.0 + pow(V/H,2))))/w + (cb/EA)*(Lu-V/w) - 1.0/w;
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return ((1.0 + (V/H)/sqrt(1.0 + pow(V/H,2)))/((V/H) + sqrt(1.0 + pow(V/H,2))))/w + (H/(w*EA)) - 1.0/w;
  };
};


double jacobian_dzdh_contact(const double V, const double H, const double w, const double Lu, const double EA, const double cb)
{
  // @todo remove this conditional statement because the equations are the same regardless of the outcome
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return (sqrt(1.0 + pow(V/H,2)) - 1.0 - (V/H)*(V/H)/sqrt(1.0 + pow(V/H,2)))/w;
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return (sqrt(1.0 + pow(V/H,2)) - 1.0 - (V/H)*(V/H)/sqrt(1.0 + pow(V/H,2)))/w;
  };
};



double jacobian_dzdv_contact(const double V, const double H, const double w, const double Lu, const double EA, const double cb)
{
  // @todo remove this conditional statement because the equations are the same regardless of the outcome
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return ((V/H)/sqrt(1.0 + pow(V/H,2)))/w + (V/(w*EA));
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return ((V/H)/sqrt(1.0 + pow(V/H,2)))/w + (V/(w*EA));
  };  
};


MAP_ERROR_CODE forward_difference_jacobian(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_ConstraintStateType_t* z_type, Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr)
{
  OuterSolveAttributes* ns = &domain->outer_loop;
  MAP_ERROR_CODE success = MAP_SAFE;
  double original_displacement = 0.0;
  const int THREE = 3;
  const int z_size = z_type->z_Len; 
  // const int m = THREE*(other_type->Fz_connect_Len); // rows
  const int n = THREE*(z_type->z_Len);              // columns
  int i = 0;
  int j = 0;
  
  for (i=0 ; i<z_size ; i++) {
    ns->b[THREE*i] = other_type->Fx_connect[i];
    ns->b[THREE*i+1] = other_type->Fy_connect[i];
    ns->b[THREE*i+2] = other_type->Fz_connect[i];      
  }

  /* First store the connect node force before applying ns->epsilon displacement */
  for (i=0 ; i<n ; i++) {            
    for (j=0 ; j<z_size ; j++) {            
      ns->jac[THREE*j][i] = -other_type->Fx_connect[j];
      ns->jac[THREE*j+1][i] = -other_type->Fy_connect[j];
      ns->jac[THREE*j+2][i] = -other_type->Fz_connect[j];
    };
  };
    
  for (j=0 ; j<z_size ; j++) {   // column               
    for (i=0 ; i<z_size ; i++) { // rows           
      original_displacement = z_type->x[j];
      z_type->x[j] += ns->epsilon;
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Forward difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j] += other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j] += other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j] += other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j] /= ns->epsilon;
      z_type->x[j] = original_displacement;
        
      original_displacement = z_type->y[j];
      z_type->y[j] += ns->epsilon;
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Forward difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+1] += other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+1] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j+1] += other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+1] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j+1] += other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+1] /= ns->epsilon;
      z_type->y[j] = original_displacement;
    
      original_displacement = z_type->z[j];
      z_type->z[j] += ns->epsilon;
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Forward difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+2] += other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+2] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j+2] += other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+2] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j+2] += other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+2] /= ns->epsilon;
      z_type->z[j] = original_displacement;
    };
  };

  /* read flag to set scaling parameter */
  if (ns->pg) {
    for (i=0 ; i<THREE*z_size ; i++) { 
      ns->jac[i][i] += (ns->ds/pow(ns->iteration_count,1.5)+ns->d);
    };
  };

  return MAP_SAFE;
};


MAP_ERROR_CODE backward_difference_jacobian(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_ConstraintStateType_t* z_type, Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr)
{
  OuterSolveAttributes* ns = &domain->outer_loop;
  MAP_ERROR_CODE success = MAP_SAFE;
  double original_displacement = 0.0;
  const int THREE = 3;
  const int z_size = z_type->z_Len; // N
  // const int m = THREE*(other_type->Fz_connect_Len); // rows
  const int n = THREE*(z_type->z_Len);              // columns
  int i = 0;
  int j = 0;
  
  for (i=0 ; i<z_size ; i++) {
    ns->b[THREE*i] = other_type->Fx_connect[i];
    ns->b[THREE*i+1] = other_type->Fy_connect[i];
    ns->b[THREE*i+2] = other_type->Fz_connect[i];      
  }

  /* First store the connect node force before applying ns->epsilon displacement */
  for (i=0 ; i<n ; i++) {            
    for (j=0 ; j<z_size ; j++) {            
      ns->jac[THREE*j][i] = other_type->Fx_connect[j];
      ns->jac[THREE*j+1][i] = other_type->Fy_connect[j];
      ns->jac[THREE*j+2][i] = other_type->Fz_connect[j];
    };
  };

  for (j=0 ; j<z_size ; j++) {   // column               
    for (i=0 ; i<z_size ; i++) { // rows           
      original_displacement = z_type->x[j];
      z_type->x[j] -= ns->epsilon;
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Backward difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j] -= other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j] -= other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j] -= other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j] /= ns->epsilon;
      z_type->x[j] = original_displacement;
        
      original_displacement = z_type->y[j];
      z_type->y[j] -= ns->epsilon;
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Backward difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+1] -= other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+1] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j+1] -= other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+1] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j+1] -= other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+1] /= ns->epsilon;
      z_type->y[j] = original_displacement;
    
      original_displacement = z_type->z[j];
      z_type->z[j] -= ns->epsilon;
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Backward difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+2] -= other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+2] /= ns->epsilon;
      ns->jac[THREE*i+1][THREE*j+2] -= other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+2] /= ns->epsilon;
      ns->jac[THREE*i+2][THREE*j+2] -= other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+2] /= ns->epsilon;
      z_type->z[j] = original_displacement;
    };
  };

  /* read flag to set scaling parameter */
  if (ns->pg) {
    for (i=0 ; i<THREE*z_size ; i++) { 
      ns->jac[i][i] += (ns->ds/pow(ns->iteration_count,1.5)+ns->d);
    };
  };

  return MAP_SAFE;
};


MAP_ERROR_CODE central_difference_jacobian(MAP_OtherStateType_t* other_type, MAP_ParameterType_t* p_type, MAP_ConstraintStateType_t* z_type, Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr)
{
  OuterSolveAttributes* ns = &domain->outer_loop;
  MAP_ERROR_CODE success = MAP_SAFE;
  double original_displacement = 0.0;
  const int THREE = 3;
  const int z_size = z_type->z_Len; //M
  // const int m = THREE*(other_type->Fz_connect_Len); // rows
  // const int n = THREE*(z_type->z_Len);              // columns
  int i = 0;
  int j = 0;
  
  for (i=0 ; i<z_size ; i++) {
    ns->b[THREE*i] = other_type->Fx_connect[i];
    ns->b[THREE*i+1] = other_type->Fy_connect[i];
    ns->b[THREE*i+2] = other_type->Fz_connect[i];      
  }

  /* First store the connect node force before applying ns->epsilon displacement */
  for (j=0 ; j<z_size ; j++) {   // column               
    for (i=0 ; i<z_size ; i++) { // rows           
      original_displacement = z_type->x[j];
      z_type->x[j] += ns->epsilon;
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Central difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j] = other_type->Fx_connect[i];      
      ns->jac[THREE*i+1][THREE*j] = other_type->Fy_connect[i];
      ns->jac[THREE*i+2][THREE*j] = other_type->Fz_connect[i];
      z_type->x[j] = original_displacement;
        
      original_displacement = z_type->y[j];
      z_type->y[j] += ns->epsilon;
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Central difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
      ns->jac[THREE*i][THREE*j+1] = other_type->Fx_connect[i];
      ns->jac[THREE*i+1][THREE*j+1] = other_type->Fy_connect[i];
      ns->jac[THREE*i+2][THREE*j+1] = other_type->Fz_connect[i];
      z_type->y[j] = original_displacement;
    
      original_displacement = z_type->z[j];
      z_type->z[j] += ns->epsilon;
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Central difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
      ns->jac[THREE*i][THREE*j+2] = other_type->Fx_connect[i];
      ns->jac[THREE*i+1][THREE*j+2] = other_type->Fy_connect[i];
      ns->jac[THREE*i+2][THREE*j+2] = other_type->Fz_connect[i];
      z_type->z[j] = original_displacement;
    };
  };
    
  for (j=0 ; j<z_size ; j++) {   // column               
    for (i=0 ; i<z_size ; i++) { // rows           
      original_displacement = z_type->x[j];
      z_type->x[j] -= ns->epsilon;
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Central difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
      ns->jac[THREE*i][THREE*j] -= other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j] /= (2*ns->epsilon);
      ns->jac[THREE*i+1][THREE*j] -= other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j] /= (2*ns->epsilon);
      ns->jac[THREE*i+2][THREE*j] -= other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j] /= (2*ns->epsilon);
      z_type->x[j] = original_displacement;
        
      original_displacement = z_type->y[j];
      z_type->y[j] -= ns->epsilon;
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Central difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+1] -= other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+1] /= (2*ns->epsilon);
      ns->jac[THREE*i+1][THREE*j+1] -= other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+1] /= (2*ns->epsilon);
      ns->jac[THREE*i+2][THREE*j+1] -= other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+1] /= (2*ns->epsilon);
      z_type->y[j] = original_displacement;
    
      original_displacement = z_type->z[j];
      z_type->z[j] -= ns->epsilon;
      success = line_solve_sequence(domain, p_type, 0.0, map_msg, ierr);
      if (success) {
        set_universal_error_with_message(map_msg, ierr, MAP_FATAL_78, "Central difference, x[%d]+delta, row %d, col %d.", j+1, THREE*i, THREE*j);
        return MAP_FATAL;
      };
      ns->jac[THREE*i][THREE*j+2] -= other_type->Fx_connect[i];
      ns->jac[THREE*i][THREE*j+2] /= (2*ns->epsilon);
      ns->jac[THREE*i+1][THREE*j+2] -= other_type->Fy_connect[i];
      ns->jac[THREE*i+1][THREE*j+2] /= (2*ns->epsilon);
      ns->jac[THREE*i+2][THREE*j+2] -= other_type->Fz_connect[i];
      ns->jac[THREE*i+2][THREE*j+2] /= (2*ns->epsilon);
      z_type->z[j] = original_displacement;
    };
  };

  /* read flag to set scaling parameter */
  if (ns->pg) {
    ns->coef = pow(ns->iteration_count,1.5);
    for (i=0 ; i<THREE*z_size ; i++) { 
      ns->jac[i][i] += (ns->ds/ns->coef + ns->d);
    };
  };

  return MAP_SAFE;
};
