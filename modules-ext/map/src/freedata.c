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


/**
 * @file 
 */


#include "freedata.h"



MAP_EXTERNCALL void MAP_InitInput_Delete(InitializationData* init_data)
{
  MAPFREE(init_data); 
};


MAP_EXTERNCALL void MAP_OtherState_Delete(Domain* domain)
{
  MAPFREE(domain);
};


MAP_ERROR_CODE free_outlist(Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr)
{
  VarTypePtr* vartype_ptr = NULL;

  if (domain->y_list) { /* if allocated, then proceed with free'ing */
    list_iterator_start(&domain->y_list->out_list_ptr);
    while (list_iterator_hasnext(&domain->y_list->out_list_ptr)) { 
      vartype_ptr = (VarTypePtr*)list_iterator_next(&domain->y_list->out_list_ptr);
      bdestroy(vartype_ptr->name);
      bdestroy(vartype_ptr->units);
    };
    list_iterator_stop(&domain->y_list->out_list_ptr);     
    
    // @rm y_list->out_list no longer exists/is useful ?
    list_destroy(&domain->y_list->out_list);     /* destroy output lists for writting information to output file */
    list_destroy(&domain->y_list->out_list_ptr); /* destroy output lists for writting information to output file */
  };
  MAPFREE(domain->y_list);
  return MAP_SAFE;
};


MAP_ERROR_CODE free_cable_library(list_t* restrict library)
{
  CableLibrary* iter_library = NULL;
  list_iterator_start(library);          /* starting an iteration "session" */
  while (list_iterator_hasnext(library)) { /* tell whether more values available */
    iter_library = (CableLibrary*)list_iterator_next(library);
    bdestroy(iter_library->label);
  };
  list_iterator_stop(library);
  return MAP_SAFE;
};


MAP_ERROR_CODE free_update_list (list_t* restrict ref_list)
{
  ReferencePoint* point_iter = NULL;
  list_iterator_start(ref_list); /* starting an iteration "session" */
  while (list_iterator_hasnext(ref_list)) { /* tell whether more values available */
    point_iter = (ReferencePoint*)list_iterator_next(ref_list);
    point_iter->x = NULL;
    point_iter->y = NULL;
    point_iter->z = NULL;
  };
  list_iterator_stop(ref_list); 

  return MAP_SAFE;
};


MAP_ERROR_CODE free_line(list_t* restrict line) 
{
  Line* line_iter = NULL;
  list_iterator_start(line); /* starting an iteration "session" */
  while (list_iterator_hasnext(line)) { /* tell whether more values available */
    line_iter = (Line*)list_iterator_next(line);
    // @rm  bdestroy(line_iter->psi.name); 
    // @rm  bdestroy(line_iter->psi.units);
    // @rm  bdestroy(line_iter->alpha.name);
    // @rm  bdestroy(line_iter->alpha.units);
    // @rm  bdestroy(line_iter->alpha_at_anchor.name);
    // @rm  bdestroy(line_iter->alpha_at_anchor.units);
    // @rm  bdestroy(line_iter->l.name); 
    // @rm  bdestroy(line_iter->l.units);
    // @rm  bdestroy(line_iter->Lb.name); 
    // @rm  bdestroy(line_iter->Lb.units);
    bdestroy(line_iter->Lu.name); 
    bdestroy(line_iter->Lu.units);
    // @rm  bdestroy(line_iter->h.name); 
    // @rm  bdestroy(line_iter->h.units);
    bdestroy(line_iter->H.name); 
    bdestroy(line_iter->H.units);
    bdestroy(line_iter->V.name); 
    bdestroy(line_iter->V.units);
    // @rm  bdestroy(line_iter->H_at_anchor.name); 
    // @rm  bdestroy(line_iter->H_at_anchor.units);
    // @rm  bdestroy(line_iter->V_at_anchor.name); 
    // @rm  bdestroy(line_iter->V_at_anchor.units);
    // @rm  bdestroy(line_iter->force_at_fairlead.fx.name); 
    // @rm  bdestroy(line_iter->force_at_fairlead.fx.units);
    // @rm  bdestroy(line_iter->force_at_fairlead.fy.name);
    // @rm  bdestroy(line_iter->force_at_fairlead.fy.units);
    // @rm  bdestroy(line_iter->force_at_fairlead.fz.name); 
    // @rm  bdestroy(line_iter->force_at_fairlead.fz.units);
    // @rm  bdestroy(line_iter->force_at_anchor.fx.name);
    // @rm  bdestroy(line_iter->force_at_anchor.fx.units);
    // @rm  bdestroy(line_iter->force_at_anchor.fy.name);
    // @rm  bdestroy(line_iter->force_at_anchor.fy.units);
    // @rm  bdestroy(line_iter->force_at_anchor.fz.name); 
    // @rm  bdestroy(line_iter->force_at_anchor.fz.units);
    // @rm  bdestroy(line_iter->T.name);
    // @rm  bdestroy(line_iter->T.units);
    // @rm  bdestroy(line_iter->tension_at_anchor.name);
    // @rm  bdestroy(line_iter->tension_at_anchor.units);
  
    /* don't let any pointers dangle */
    line_iter->line_property = NULL;      
    line_iter->label = NULL;
    line_iter->line_tension = NULL;
    line_iter->anchor = NULL; 
    line_iter->fairlead = NULL;
  };
  list_iterator_stop(line); /* ending the iteration "session" */  
  return MAP_SAFE;
};


MAP_ERROR_CODE free_node(list_t *restrict node)
{
  Node* iterNode = NULL;
  MAP_ERROR_CODE success = MAP_SAFE;
  list_iterator_start(node);            /* starting an iteration "session" */
  while (list_iterator_hasnext(node)) { /* tell whether more values available */ 
    iterNode = (Node*)list_iterator_next(node);

    success = bdestroy(iterNode->M_applied.name); 
    success = bdestroy(iterNode->M_applied.units);
    success = bdestroy(iterNode->B_applied.name);
    success = bdestroy(iterNode->B_applied.units);

    success = bdestroy(iterNode->external_force.fx.name);
    success = bdestroy(iterNode->external_force.fx.units);
    success = bdestroy(iterNode->external_force.fy.name);
    success = bdestroy(iterNode->external_force.fy.units);
    success = bdestroy(iterNode->external_force.fz.name);
    success = bdestroy(iterNode->external_force.fz.units);

    success = bdestroy(iterNode->position_ptr.x.name); 
    success = bdestroy(iterNode->position_ptr.x.units);
    success = bdestroy(iterNode->position_ptr.y.name); 
    success = bdestroy(iterNode->position_ptr.y.units);
    success = bdestroy(iterNode->position_ptr.z.name); 
    success = bdestroy(iterNode->position_ptr.z.units);

    success = bdestroy(iterNode->sum_force_ptr.fx.name); 
    success = bdestroy(iterNode->sum_force_ptr.fx.units);
    success = bdestroy(iterNode->sum_force_ptr.fy.name); 
    success = bdestroy(iterNode->sum_force_ptr.fy.units);
    success = bdestroy(iterNode->sum_force_ptr.fz.name); 
    success = bdestroy(iterNode->sum_force_ptr.fz.units);
  };
  list_iterator_stop(node); /* ending the iteration "session" */  
  return MAP_SAFE;
};


MAP_ERROR_CODE free_vessel(Vessel* floater) 
{
  /* Now delete the vessel information */
  MAPFREE(floater->xi);
  MAPFREE(floater->yi);
  MAPFREE(floater->zi);

  bdestroy(floater->displacement.x.name);
  bdestroy(floater->displacement.x.units);
  bdestroy(floater->displacement.y.name);
  bdestroy(floater->displacement.y.units);
  bdestroy(floater->displacement.z.name);
  bdestroy(floater->displacement.z.units);
  
  bdestroy(floater->ref_origin.x.name);
  bdestroy(floater->ref_origin.x.units);
  bdestroy(floater->ref_origin.y.name);
  bdestroy(floater->ref_origin.y.units);
  bdestroy(floater->ref_origin.z.name);
  bdestroy(floater->ref_origin.z.units);
          
  bdestroy(floater->line_sum_force.fx.name);
  bdestroy(floater->line_sum_force.fx.units);
  bdestroy(floater->line_sum_force.fy.name);
  bdestroy(floater->line_sum_force.fy.units);
  bdestroy(floater->line_sum_force.fz.name);
  bdestroy(floater->line_sum_force.fz.units);

  bdestroy(floater->orientation.phi.name);
  bdestroy(floater->orientation.phi.units);
  bdestroy(floater->orientation.the.name);
  bdestroy(floater->orientation.the.units);
  bdestroy(floater->orientation.psi.name);
  bdestroy(floater->orientation.psi.units);
  return MAP_SAFE;
}


MAP_ERROR_CODE map_free_types(MAP_InputType_t* u_type, MAP_ParameterType_t* p_type, MAP_ContinuousStateType_t* x_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type)
{
  /* inputs */
  MAPFREE(u_type->x);
  MAPFREE(u_type->y);
  MAPFREE(u_type->z);

  /* parameters are skipped for now; they are set in fortran since depth, gravity and sea density are set by glue code */

  /* continuous state */

  /* constraint state */  
  MAPFREE(z_type->H);     
  MAPFREE(z_type->V);     
  MAPFREE(z_type->x);     
  MAPFREE(z_type->y);     
  MAPFREE(z_type->z);     

  /* other state */
  MAPFREE(other_type->H); 
  MAPFREE(other_type->V); 
  MAPFREE(other_type->Ha);
  MAPFREE(other_type->Va);
  MAPFREE(other_type->x); 
  MAPFREE(other_type->y); 
  MAPFREE(other_type->z); 
  MAPFREE(other_type->xa);
  MAPFREE(other_type->ya);
  MAPFREE(other_type->za);
  MAPFREE(other_type->Fx_connect); 
  MAPFREE(other_type->Fy_connect); 
  MAPFREE(other_type->Fz_connect); 
  MAPFREE(other_type->Fx_anchor); 
  MAPFREE(other_type->Fy_anchor); 
  MAPFREE(other_type->Fz_anchor); 

  /* outputs */
  MAPFREE(y_type->Fx);    
  MAPFREE(y_type->Fy);    
  MAPFREE(y_type->Fz);    
  MAPFREE(y_type->wrtOutput);
  MAPFREE(y_type->WriteOutput);
  
  return MAP_SAFE;
};


MAP_ERROR_CODE free_outer_solve_data(OuterSolveAttributes* ns, const int size, char* map_msg, MAP_ERROR_CODE* ierr)
{
  // const int N = ns->max_krylov_its + 1;
  const int SIZE = 3*size;
  int i = 0;

  if (ns->jac) { /* is it allocated? */
    for(i=0 ; i<SIZE ; i++) {
      MAPFREE(ns->jac[i]);
    };
  };
 
  if (ns->l) { /* is it allocated? */
    for(i=0 ; i<SIZE ; i++) {
     MAPFREE(ns->l[i]);
    };
  };

  if (ns->u) { /* is it allocated? */
    for(i=0 ; i<SIZE ; i++) {
      MAPFREE(ns->u[i]);
   };  
  };

  if (ns->V) { /* is it allocated? */
    for(i=0 ; i<SIZE ; i++) {
      MAPFREE(ns->V[i]);
   };  
  };

  if (ns->AV) { /* is it allocated? */
    for(i=0 ; i<SIZE ; i++) {
      MAPFREE(ns->AV[i]);
   };  
  };


  MAPFREE(ns->jac);
  MAPFREE(ns->AV);
  MAPFREE(ns->av);
  MAPFREE(ns->V);
  MAPFREE(ns->l);
  MAPFREE(ns->u);
  MAPFREE(ns->b);
  MAPFREE(ns->w);
  MAPFREE(ns->q);
  MAPFREE(ns->x);  
  MAPFREE(ns->y);
  MAPFREE(ns->C);
  return MAP_SAFE;
};

